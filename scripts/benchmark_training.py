#!/usr/bin/env python3
"""Time real ESM3 adapter optimizer steps on prepared complexes; discard model updates."""

import argparse
import json
import os
import platform
import statistics
import subprocess
import time
from pathlib import Path

import numpy as np
import psutil
import torch

from mertk.dataset import load_batch, read_records
from mertk.runtime import resolve_device, resolve_precision, synchronize
from mertk.training import backward_microbatch, build_model


def system_info(device):
    memory = psutil.virtual_memory()
    info = {
        "platform": platform.platform(),
        "torch": str(torch.__version__),
        "cpu_cores": psutil.cpu_count(),
        "ram_gib": memory.total / 2**30,
        "available_ram_gib": memory.available / 2**30,
        "mps_cpu_fallback_env": os.environ.get("PYTORCH_ENABLE_MPS_FALLBACK", "0"),
    }
    try:
        info["swap_used_gib"] = psutil.swap_memory().used / 2**30
    except OSError:
        pass
    if platform.system() == "Darwin":
        result = subprocess.run(
            ["sysctl", "-n", "machdep.cpu.brand_string"], capture_output=True, text=True
        )
        if result.returncode == 0:
            info["chip"] = result.stdout.strip()
    if device.type == "mps":
        info["mps_recommended_max_gib"] = torch.mps.recommended_max_memory() / 2**30
    return info


def memory_snapshot(device):
    result = {
        "process_rss_gib": psutil.Process().memory_info().rss / 2**30,
        "available_ram_gib": psutil.virtual_memory().available / 2**30,
    }
    if device.type == "mps":
        result.update(
            mps_tensors_gib=torch.mps.current_allocated_memory() / 2**30,
            mps_driver_gib=torch.mps.driver_allocated_memory() / 2**30,
        )
    elif device.type == "cuda":
        result.update(cuda_peak_gib=torch.cuda.max_memory_allocated(device) / 2**30)
    return result


def representative_paths(manifest):
    by_length = {}
    for row in read_records(manifest):
        with np.load(row["path"], allow_pickle=False) as example:
            length = len(str(example["binder_sequence"])) + len(str(example["target_sequence"])) + 3
        by_length.setdefault(length, row["path"])
    lengths = sorted(by_length)
    return [by_length[length] for length in sorted({lengths[0], lengths[-1]})]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("configs/train.json"))
    parser.add_argument("--device", choices=["auto", "mps", "cpu", "cuda"], default="auto")
    parser.add_argument("--precision", choices=["auto", "fp32", "bf16"], default="auto")
    parser.add_argument("--examples", nargs="+", type=Path)
    parser.add_argument("--steps", type=int, default=3, help="Measured optimizer steps per example")
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--accumulation-steps", type=int, default=1)
    parser.add_argument("--cpu-threads", type=int, default=6)
    parser.add_argument(
        "--gradient-checkpointing", action=argparse.BooleanOptionalAction, default=True
    )
    parser.add_argument("--objective", choices=["masked", "denoise"], default="masked")
    parser.add_argument("--output", type=Path, default=Path("outputs/training_benchmark.json"))
    args = parser.parse_args()
    if min(args.steps, args.accumulation_steps, args.cpu_threads) < 1 or args.warmup < 1:
        parser.error("steps, accumulation, CPU threads, and warmup must be positive")
    if args.output.exists():
        parser.error(f"Output already exists: {args.output}; choose another path")
    device = resolve_device(args.device)
    precision = resolve_precision(device, args.precision)
    config = json.loads(args.config.read_text())
    config.update(
        device=str(device),
        precision=precision,
        objective=args.objective,
        gradient_checkpointing=args.gradient_checkpointing,
    )
    paths = args.examples or representative_paths(config["manifest"])
    torch.set_num_threads(args.cpu_threads)
    torch.manual_seed(config["seed"])
    mask_rng = torch.Generator().manual_seed(config["seed"])
    result = {
        "system": system_info(device),
        "device": str(device),
        "precision": precision,
        "gradient_checkpointing": args.gradient_checkpointing,
        "objective": args.objective,
        "accumulation_steps": args.accumulation_steps,
        "warmup_steps": args.warmup,
        "cpu_threads": args.cpu_threads,
        "performance_only": True,
        "checkpoint_saved": False,
        "examples": [],
    }
    print(json.dumps(result, indent=2), flush=True)
    started = time.perf_counter()
    model = build_model(config, device)
    synchronize(device)
    result["load_seconds"] = time.perf_counter() - started
    parameters = [p for p in model.parameters() if p.requires_grad]
    model.train()
    for path in paths:
        batch = load_batch(path, model.tokenizer, device, require_structure=True)
        optimizer = torch.optim.AdamW(parameters, lr=config["learning_rate"], weight_decay=0.01)
        timings, losses, memory = [], [], []
        if device.type == "cuda":
            torch.cuda.reset_peak_memory_stats(device)
        for step in range(args.warmup + args.steps):
            synchronize(device)
            started = time.perf_counter()
            optimizer.zero_grad(set_to_none=True)
            loss = 0.0
            for _ in range(args.accumulation_steps):
                loss += (
                    backward_microbatch(
                        model, batch, device, 0.5, mask_rng, args.accumulation_steps, precision
                    )
                    / args.accumulation_steps
                )
            torch.nn.utils.clip_grad_norm_(parameters, 1.0, error_if_nonfinite=True)
            optimizer.step()
            synchronize(device)
            elapsed = time.perf_counter() - started
            measured = step >= args.warmup
            print(
                json.dumps(
                    {
                        "example": str(path),
                        "tokens": batch["sequence_tokens"].shape[1],
                        "step": step + 1,
                        "measured": measured,
                        "seconds": elapsed,
                        "loss": loss,
                    }
                ),
                flush=True,
            )
            if measured:
                timings.append(elapsed)
                losses.append(loss)
                memory.append(memory_snapshot(device))
        if any(p.grad is not None for p in model.parameters() if not p.requires_grad):
            raise RuntimeError("A frozen parameter unexpectedly received a gradient")
        changed = any(
            p.detach().abs().sum().item() > 0
            for name, p in model.named_parameters()
            if name.endswith("lora_b")
        )
        if config["use_lora"] and not changed:
            raise RuntimeError("LoRA matrices did not update")
        median = statistics.median(timings)
        result["examples"].append(
            {
                "path": str(path),
                "tokens": batch["sequence_tokens"].shape[1],
                "seconds_per_optimizer_step": timings,
                "median_seconds": median,
                "complexes_per_second": args.accumulation_steps / median,
                "estimated_1000_steps_hours_excluding_eval": median * 1000 / 3600,
                "losses": losses,
                "lora_update_verified": changed,
                "memory_after_step": memory,
            }
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"Benchmark saved to {args.output}. Model updates discarded.", flush=True)


if __name__ == "__main__":
    main()
