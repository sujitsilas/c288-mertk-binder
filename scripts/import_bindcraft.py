#!/usr/bin/env python3
import argparse
from pathlib import Path

from mertk.bindcraft import import_complexes

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Import accepted BindCraft complexes and retain original-sequence controls"
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="An accepted complex PDB or the Accepted directory",
    )
    parser.add_argument("--output", type=Path, default=Path("data/application"))
    parser.add_argument("--binder-chain", default="B")
    parser.add_argument("--target-chain", default="A")
    args = parser.parse_args()
    paths = sorted(args.input.glob("*.pdb")) if args.input.is_dir() else [args.input]
    rows = import_complexes(paths, args.output, args.binder_chain, args.target_chain)
    print(
        f"Imported {len(rows)} application holdouts; original sequences saved as bindcraft_original.fasta"
    )
