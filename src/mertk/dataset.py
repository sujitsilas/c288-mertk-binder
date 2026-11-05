"""Single-complex batches retain full chain sequences and explicit chain boundaries."""

import json
from pathlib import Path

import numpy as np
import torch


def read_records(path, split=None):
    rows = [json.loads(line) for line in Path(path).read_text().splitlines() if line.strip()]
    for row in rows:
        if "path" in row and not Path(row["path"]).is_absolute():
            row["path"] = str((Path(path).resolve().parent / row["path"]).resolve())
    if split is not None:
        rows = [row for row in rows if row["split"] == split]
    if not rows:
        raise ValueError(f"No examples for split {split!r} in {path}")
    return rows


def make_batch(example, tokenizer, device="cpu", require_structure=False):
    binder, target = str(example["binder_sequence"]), str(example["target_sequence"])
    tokens = torch.tensor(tokenizer.encode(binder + "|" + target), dtype=torch.long)
    n, m, k = len(binder), len(target), example["neighbor_xyz"].shape[1]
    length = n + m + 3
    if len(tokens) != length:
        raise ValueError("Sequence tokenization is not one token per residue")
    backbone = torch.full((length, 3, 3), float("nan"))
    backbone[1 : n + 1] = torch.from_numpy(example["binder_backbone"].copy())
    backbone[n + 2 : -1] = torch.from_numpy(example["target_backbone"].copy())
    design_mask = torch.zeros(length, dtype=torch.bool)
    design_mask[1 : n + 1] = torch.from_numpy(example["design_mask"].copy())
    if not design_mask.any():
        raise ValueError("Example has no designable interface residues")
    if not torch.isfinite(backbone[design_mask]).all():
        raise ValueError("Designable residues require finite backbone coordinates")
    batch = {
        "sequence_tokens": tokens,
        "backbone": backbone,
        "design_mask": design_mask,
        "chain_id": torch.tensor([0] * (n + 2) + [1] * (m + 1)),
        "neighbor_xyz": torch.zeros(length, k, 3),
        "neighbor_features": torch.zeros(length, k, 3, dtype=torch.long),
        "neighbor_mask": torch.zeros(length, k, dtype=torch.bool),
    }
    for name in ("neighbor_xyz", "neighbor_features", "neighbor_mask"):
        batch[name][1 : n + 1] = torch.from_numpy(example[name].copy())
    if "binder_structure_tokens" in example and "target_structure_tokens" in example:
        structure = torch.full((length,), 4096, dtype=torch.long)
        structure[1 : n + 1] = torch.from_numpy(example["binder_structure_tokens"].copy())
        structure[n + 2 : -1] = torch.from_numpy(example["target_structure_tokens"].copy())
        batch["structure_tokens"] = structure
    elif require_structure:
        raise ValueError(
            "Cache backbone structure tokens with scripts/cache_structure_tokens.py first"
        )
    return {name: value.unsqueeze(0).to(device) for name, value in batch.items()}


def load_batch(path, tokenizer, device="cpu", require_structure=False):
    with np.load(path, allow_pickle=False) as example:
        return make_batch(example, tokenizer, device, require_structure)


def synthetic_example():
    """Non-biological fixture for testing tensor plumbing, never training evidence."""
    n, m, k = 12, 10, 4
    ca = np.stack((np.arange(n) * 3.8, np.zeros(n), np.zeros(n)), -1).astype(np.float32)
    backbone = np.stack((ca + [-1, 1, 0], ca, ca + [1, 0, 0]), 1).astype(np.float32)
    return {
        "binder_sequence": np.array("ACDEFGHIKLMN"),
        "target_sequence": np.array("PQRSTVWYAC"),
        "binder_backbone": backbone,
        "target_backbone": backbone[:m] + [0, 5, 0],
        "neighbor_xyz": np.repeat((ca + [0, 4, 0])[:, None, :], k, axis=1).astype(np.float32),
        "neighbor_features": np.zeros((n, k, 3), dtype=np.int64),
        "neighbor_mask": np.ones((n, k), dtype=bool),
        "design_mask": np.ones(n, dtype=bool),
        "binder_oxygen": (backbone[:, 2] + [0.3, 1.0, 0]).astype(np.float32),
        "target_oxygen": (backbone[:m, 2] + [0.3, 6.0, 0]).astype(np.float32),
    }
