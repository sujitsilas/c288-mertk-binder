"""Single-GPU training with adapter-only checkpoints and deterministic evaluation."""

import hashlib
import json
import random
import subprocess
import time
from pathlib import Path

import torch

from mertk.dataset import load_batch, make_batch, read_records, synthetic_example
from mertk.model import InterfaceModel, load_backbone, weight_identity
from mertk.objectives import corrupt, reconstruction_loss
from mertk.runtime import (
    accelerator_rng_state,
    autocast,
    resolve_device,
    resolve_precision,
    restore_accelerator_rng,
    synchronize,
)


def build_model(config, device, tiny=False):
    if config.get("model_type", "esm3") == "atom_decoder":
        from mertk.atom_decoder import AtomDecoder

        return AtomDecoder(config["hidden"]).to(device)
    if config.get("model_type", "esm3") != "esm3":
        raise ValueError("Unknown model_type")
    base = load_backbone(config["weights"], device, tiny=tiny)
    return InterfaceModel(
        base,
        use_lora=config["use_lora"],
        use_geometry=config["use_geometry"],
        use_time=config["objective"] == "denoise",
        rank=config["rank"],
        alpha=config["alpha"],
        dropout=config["dropout"],
        hidden=config["hidden"],
        gradient_checkpointing=config["gradient_checkpointing"],
        atom_context=config.get("atom_context", "all"),
    ).to(device)


def get_batch(row, model, device, tiny=False):
    tokenizer = model.tokenizer
    if tiny:
        return make_batch(synthetic_example(), tokenizer, device)
    return load_batch(row["path"], tokenizer, device, require_structure=True)


def audit_splits(rows):
    """Refuse unassigned data or any cross-split component/entry/sequence overlap."""
    ownership = {}
    for row in rows:
        split = row["split"]
        if split in {"unsplit", ""}:
            raise ValueError("Cluster and split examples before training")
        if split == "excluded_homologous_pair":
            continue
        if split not in {"train", "val", "test", "application_holdout"}:
            raise ValueError(f"Unknown split: {split}")
        for name in ("component", "pdb_id", "binder_key", "target_key"):
            if name not in row:
                raise ValueError(f"Missing split audit field: {name}")
            key = ("sequence" if name.endswith("_key") else name, row[name])
            if key in ownership and ownership[key] != split:
                raise ValueError(f"Data leakage across splits: {key}")
            ownership[key] = split


def backward_microbatch(
    model, batch, device, probability, mask_rng, accumulation_steps=1, precision="auto"
):
    tokens, selected = corrupt(
        batch["sequence_tokens"],
        batch["design_mask"],
        model.tokenizer.mask_token_id,
        probability,
        mask_rng,
    )
    with autocast(device, precision):
        logits = model(
            batch, sequence_tokens=tokens, time=torch.tensor([probability], device=device)
        )
        loss = reconstruction_loss(logits, batch["sequence_tokens"], selected)
    if not torch.isfinite(loss):
        raise FloatingPointError("Non-finite training loss")
    (loss / accumulation_steps).backward()
    return loss.item()


@torch.no_grad()
def evaluate(model, rows, device, tiny=False, seed=2026, precision="auto"):
    model.eval()
    total_loss, correct, total = 0.0, 0, 0
    tokenizer = model.tokenizer
    from mertk.data import AA

    aa_ids = torch.tensor([tokenizer.convert_tokens_to_ids(a) for a in AA], device=device)
    for row in rows:
        batch = get_batch(row, model, device, tiny)
        # Full-interface reconstruction is identical for every ablation.
        tokens, selected = corrupt(
            batch["sequence_tokens"],
            batch["design_mask"],
            tokenizer.mask_token_id,
            probability=1.0,
            generator=torch.Generator().manual_seed(seed),
        )
        with autocast(device, precision):
            logits = model(batch, sequence_tokens=tokens, time=torch.ones(1, device=device))
        count = int(selected.sum())
        loss = reconstruction_loss(logits, batch["sequence_tokens"], selected)
        if not 0 < count <= selected.numel() or not torch.isfinite(loss):
            raise FloatingPointError("Invalid validation mask or non-finite validation loss")
        total_loss += loss.item() * count
        predictions = aa_ids[logits[selected][:, aa_ids].argmax(-1)]
        matches = int((predictions == batch["sequence_tokens"][selected]).sum())
        if not 0 <= matches <= count:
            raise FloatingPointError("Invalid validation recovery count")
        correct += matches
        total += count
    if not total:
        raise ValueError("Evaluation has no masked residues")
    return {
        "loss": total_loss / total,
        "sequence_recovery": correct / total,
        "residues": total,
        "complexes": len(rows),
    }


def file_hash(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def train(config, smoke=False, resume=None):
    if config["objective"] not in {"masked", "denoise"}:
        raise ValueError("objective must be masked or denoise")
    if config["steps"] < 1 or config["accumulation_steps"] < 1 or config["eval_every"] < 1:
        raise ValueError("steps, accumulation_steps, and eval_every must be positive")
    device = resolve_device(
        config.get("device", "cpu") if not smoke or config.get("smoke_device") else "cpu"
    )
    precision = resolve_precision(device, config.get("precision", "auto"))
    if config.get("cpu_threads") is not None:
        if config["cpu_threads"] < 1:
            raise ValueError("cpu_threads must be positive")
        torch.set_num_threads(config["cpu_threads"])
    torch.manual_seed(config["seed"])
    rng = random.Random(config["seed"])
    mask_rng = torch.Generator().manual_seed(config["seed"])
    if smoke:
        train_rows = val_rows = [{"id": "synthetic-smoke-only"}]
    else:
        all_rows = read_records(config["manifest"])
        audit_splits(all_rows)
        train_rows = [row for row in all_rows if row["split"] == "train"]
        val_rows = [row for row in all_rows if row["split"] == "val"]
        if not train_rows or not val_rows:
            raise ValueError(
                "Nonempty train and val families are required; download more complexes"
            )
    model = build_model(config, device, tiny=smoke)
    parameters = [p for p in model.parameters() if p.requires_grad]
    if not parameters:
        raise ValueError("No trainable parameters; evaluate the frozen baseline instead")
    optimizer = torch.optim.AdamW(parameters, lr=config["learning_rate"], weight_decay=0.01)
    output = Path(config["output"])
    output.mkdir(parents=True, exist_ok=True)
    identity = (
        {"kind": "random-tiny-smoke"}
        if smoke
        else (
            {"kind": "atom_decoder"}
            if config.get("model_type") == "atom_decoder"
            else weight_identity(config["weights"])
        )
    )
    metadata = {
        "config": config,
        "smoke_test": smoke,
        "weights": identity,
        "manifest_sha256": None if smoke else file_hash(config["manifest"]),
        "trainable_parameters": sum(p.numel() for p in parameters),
        "total_parameters": sum(p.numel() for p in model.parameters()),
        "runtime": {"device": str(device), "precision": precision, "torch": str(torch.__version__)},
    }
    revision = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True)
    metadata["git_revision"] = revision.stdout.strip() if revision.returncode == 0 else None
    start = 0
    if resume:
        saved = torch.load(resume, map_location="cpu", weights_only=True)
        for key in ("config", "smoke_test", "weights", "manifest_sha256"):
            if saved["metadata"][key] != metadata[key]:
                raise ValueError(f"Resume metadata mismatch: {key}")
        if saved["metadata"].get("runtime", metadata["runtime"]) != metadata["runtime"]:
            raise ValueError("Resume runtime differs from the saved device/precision")
        model.load_adapter_state_dict(saved["adapter"])
        optimizer.load_state_dict(saved["optimizer"])
        rng.setstate(saved["python_rng"])
        mask_rng.set_state(saved["mask_rng"])
        torch.set_rng_state(saved["torch_rng"])
        restore_accelerator_rng(device, saved.get("accelerator_rng", saved.get("cuda_rng")))
        start = saved["step"]
        # Rolling back to an older saved step must not leave duplicate/future metrics.
        metrics_path = output / "metrics.jsonl"
        if metrics_path.exists():
            lines = [
                line
                for line in metrics_path.read_text().splitlines()
                if line.strip() and json.loads(line)["step"] <= start
            ]
            metrics_path.write_text("".join(line + "\n" for line in lines))
    elif (output / "last.pt").exists() or (output / "metadata.json").exists():
        raise FileExistsError(f"Existing run at {output}; use --resume or choose another output")
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(
        json.dumps(
            {
                key: metadata[key]
                for key in ["smoke_test", "trainable_parameters", "total_parameters"]
            }
        ),
        flush=True,
    )
    for step in range(start, config["steps"]):
        synchronize(device)
        step_started = time.perf_counter()
        model.train()
        optimizer.zero_grad(set_to_none=True)
        train_loss = 0.0
        lr = config["learning_rate"] * min(1.0, (step + 1) / max(1, config["warmup_steps"]))
        for group in optimizer.param_groups:
            group["lr"] = lr
        for _ in range(config["accumulation_steps"]):
            row = rng.choice(train_rows)
            batch = get_batch(row, model, device, smoke)
            probability = (
                rng.uniform(0.05, 1.0)
                if config["objective"] == "denoise"
                else config["mask_probability"]
            )
            loss = backward_microbatch(
                model, batch, device, probability, mask_rng, config["accumulation_steps"], precision
            )
            train_loss += loss / config["accumulation_steps"]
        torch.nn.utils.clip_grad_norm_(parameters, 1.0, error_if_nonfinite=True)
        optimizer.step()
        synchronize(device)
        metrics = {
            "step": step + 1,
            "train_loss": train_loss,
            "learning_rate": lr,
            "step_seconds": time.perf_counter() - step_started,
        }
        if (step + 1) % config["eval_every"] == 0 or step + 1 == config["steps"]:
            metrics["validation"] = evaluate(model, val_rows, device, smoke, precision=precision)
            saved = {
                "step": step + 1,
                "adapter": model.adapter_state_dict(),
                "optimizer": optimizer.state_dict(),
                "metadata": metadata,
                "python_rng": rng.getstate(),
                "mask_rng": mask_rng.get_state(),
                "torch_rng": torch.get_rng_state(),
                "cuda_rng": torch.cuda.get_rng_state_all() if device.type == "cuda" else [],
                "accelerator_rng": accelerator_rng_state(device),
            }
            torch.save(saved, output / "last.tmp.pt")
            (output / "last.tmp.pt").replace(output / "last.pt")
        with (output / "metrics.jsonl").open("a") as handle:
            handle.write(json.dumps(metrics) + "\n")
        print(json.dumps(metrics), flush=True)
    return output / "last.pt"
