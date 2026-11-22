#!/usr/bin/env python3
"""Import completed BindCraft, RFD4-Proteina or refined protein complexes."""

import argparse
from pathlib import Path

from mertk.candidates import import_complexes


def main(default_arm=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="PDB/mmCIF file or directory")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--arm", default=default_arm, required=default_arm is None)
    parser.add_argument("--campaign-id", required=True)
    parser.add_argument("--epitope-id", required=True)
    parser.add_argument("--target-reference", type=Path, required=True)
    parser.add_argument("--run-record", type=Path, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--binder-chain", default="B")
    parser.add_argument("--target-chain", default="A")
    args = parser.parse_args()
    paths = (
        sorted(p for p in args.input.iterdir() if p.suffix in {".pdb", ".cif"})
        if args.input.is_dir()
        else [args.input]
    )
    rows = import_complexes(
        paths,
        args.output,
        arm=args.arm,
        campaign_id=args.campaign_id,
        epitope_id=args.epitope_id,
        target_reference=args.target_reference,
        run_record=args.run_record,
        seed=args.seed,
        binder_chain=args.binder_chain,
        target_chain=args.target_chain,
    )
    print(f"Imported {len(rows)} candidates into {args.output / 'candidates.jsonl'}")


if __name__ == "__main__":
    main()
