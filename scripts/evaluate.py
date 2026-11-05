#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

import torch

from mertk.dataset import read_records
from mertk.model import weight_identity
from mertk.runtime import resolve_device
from mertk.training import audit_splits, build_model, evaluate, file_hash

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Evaluate matched full-interface recovery masks")
    parser.add_argument("--config", type=Path, default=Path("configs/train.json"))
    parser.add_argument("--checkpoint", type=Path, help="Omit for frozen ESM3 baseline")
    parser.add_argument("--split", choices=["val", "test", "application_holdout"], default="test")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--device", choices=["auto", "cpu", "mps", "cuda"])
    parser.add_argument("--precision", choices=["auto", "fp32", "bf16"])
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    for key in ("device", "precision"):
        if getattr(args, key) is not None:
            config[key] = getattr(args, key)
    saved = None
    if args.checkpoint:
        saved = torch.load(args.checkpoint, map_location="cpu", weights_only=True)
        if saved["metadata"]["smoke_test"]:
            raise ValueError(
                "A synthetic smoke checkpoint cannot be used for scientific evaluation"
            )
        for key in ["use_lora", "use_geometry", "objective", "rank", "alpha", "dropout", "hidden"]:
            config[key] = saved["metadata"]["config"][key]
        if saved["metadata"]["weights"] != weight_identity(config["weights"]):
            raise ValueError("Checkpoint base weights differ from evaluation weights")
        if saved["metadata"]["manifest_sha256"] != file_hash(config["manifest"]):
            raise ValueError("Evaluation must use the same split manifest as training")
    else:
        config.update(use_lora=False, use_geometry=False, objective="masked")
    rows = read_records(config["manifest"])
    audit_splits(rows)
    selected = [row for row in rows if row["split"] == args.split]
    device = resolve_device(config["device"])
    model = build_model(config, device)
    if saved:
        model.load_adapter_state_dict(saved["adapter"])
    metrics = evaluate(model, selected, device, precision=config.get("precision", "auto"))
    result = {
        "split": args.split,
        "checkpoint": str(args.checkpoint) if saved else "frozen_esm3",
        "config": config,
        "metrics": metrics,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(metrics, indent=2))
