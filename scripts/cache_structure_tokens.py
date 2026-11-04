#!/usr/bin/env python3
"""Encode only N/CA/C coordinates with ESM3's frozen structure encoder."""

import argparse
from pathlib import Path

import numpy as np
import torch

from mertk.dataset import read_records
from mertk.model import load_structure_encoder


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--examples", default="data/processed/examples.jsonl")
    parser.add_argument("--weights", default="weights/esm3")
    parser.add_argument("--device", default="cpu")
    args = parser.parse_args()
    encoder = load_structure_encoder(args.weights, args.device)
    for row in read_records(args.examples):
        path = Path(row["path"])
        with np.load(path, allow_pickle=False) as data:
            example = dict(data)
        with torch.no_grad():
            for chain in ("binder", "target"):
                coords = torch.tensor(example[f"{chain}_backbone"], device=args.device)[None]
                valid = torch.isfinite(coords).all(dim=(-1, -2))
                _, tokens = encoder.encode(coords)
                example[f"{chain}_structure_tokens"] = (
                    tokens.masked_fill(~valid, 4096)[0].cpu().numpy()
                )
        temporary = path.with_suffix(".tmp.npz")
        np.savez_compressed(temporary, **example)
        temporary.replace(path)
        print(f"Encoded {row['id']}", flush=True)


if __name__ == "__main__":
    main()
