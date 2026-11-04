#!/usr/bin/env python3
"""Redesign a prepared interface while preserving its target and fixed residues."""

import argparse
import json
from pathlib import Path

import numpy as np
import torch

from mertk.dataset import load_batch
from mertk.model import weight_identity
from mertk.objectives import sample_sequence
from mertk.training import autocast, build_model

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--example", type=Path, required=True, help="Prepared interface .npz")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--weights", default="weights/esm3")
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--samples", type=int, default=8)
    parser.add_argument("--steps", type=int, default=16)
    parser.add_argument("--temperature", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    if args.samples < 1:
        raise ValueError("samples must be positive")
    saved = torch.load(args.checkpoint, map_location="cpu", weights_only=True)
    if saved["metadata"]["smoke_test"]:
        raise ValueError("Use a trained real ESM3 checkpoint for design")
    config = dict(saved["metadata"]["config"], weights=args.weights)
    if saved["metadata"]["weights"] != weight_identity(args.weights):
        raise ValueError("Base weights differ from training weights")
    device = torch.device(args.device)
    model = build_model(config, device)
    model.load_adapter_state_dict(saved["adapter"])
    tokenizer = model.backbone.tokenizers.sequence
    batch = load_batch(args.example, tokenizer, device, require_structure=True)
    with np.load(args.example, allow_pickle=False) as example:
        binder_length = len(str(example["binder_sequence"]))
    sequences = []
    for index in range(args.samples):
        with autocast(device):
            tokens = sample_sequence(
                model, batch, tokenizer, args.steps, args.temperature, args.seed + index
            )
        fixed = ~batch["design_mask"]
        if not torch.equal(tokens[fixed], batch["sequence_tokens"][fixed]):
            raise RuntimeError("A fixed residue was changed")
        sequences.append(
            "".join(tokenizer.convert_ids_to_tokens(tokens[0, 1 : binder_length + 1].tolist()))
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        "".join(f">design_{i} seed={args.seed + i}\n{seq}\n" for i, seq in enumerate(sequences))
    )
    args.output.with_suffix(".json").write_text(
        json.dumps(
            {
                "checkpoint": str(args.checkpoint),
                "example": str(args.example),
                "seed": args.seed,
                "samples": args.samples,
                "unique_sequences": len(set(sequences)),
                "steps": args.steps,
                "temperature": args.temperature,
                "objective": config["objective"],
                "scored_for_binding": False,
            },
            indent=2,
        )
        + "\n"
    )
