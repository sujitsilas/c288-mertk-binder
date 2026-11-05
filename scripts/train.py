#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

from mertk.training import train

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train ESM3 LoRA / geometric interface adapters")
    parser.add_argument("--config", type=Path, default=Path("configs/train.json"))
    parser.add_argument("--output")
    parser.add_argument("--steps", type=int)
    parser.add_argument("--device", choices=["auto", "cpu", "mps", "cuda"])
    parser.add_argument("--precision", choices=["auto", "fp32", "bf16"])
    parser.add_argument("--accumulation-steps", type=int)
    parser.add_argument("--cpu-threads", type=int)
    parser.add_argument(
        "--gradient-checkpointing", action=argparse.BooleanOptionalAction, default=None
    )
    parser.add_argument("--objective", choices=["masked", "denoise"])
    parser.add_argument("--ablation", choices=["combined", "lora", "geometry"], default="combined")
    parser.add_argument(
        "--smoke-test",
        action="store_true",
        help="Tiny random ESM3 and synthetic data (CPU unless --device is given)",
    )
    parser.add_argument("--resume", type=Path)
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    for key in (
        "device",
        "precision",
        "accumulation_steps",
        "cpu_threads",
        "gradient_checkpointing",
    ):
        if getattr(args, key) is not None:
            config[key] = getattr(args, key)
    if args.smoke_test and args.device is not None:
        config["smoke_device"] = True
    if args.output:
        config["output"] = args.output
    if args.steps is not None:
        config["steps"] = args.steps
    if args.objective:
        config["objective"] = args.objective
    config["use_lora"] = args.ablation in {"combined", "lora"}
    config["use_geometry"] = args.ablation in {"combined", "geometry"}
    if args.smoke_test:
        config.update(
            steps=args.steps or 2,
            accumulation_steps=1,
            eval_every=1,
            output=args.output or "outputs/smoke",
        )
    print(f"Checkpoint: {train(config, smoke=args.smoke_test, resume=args.resume)}")
