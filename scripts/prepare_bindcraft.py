#!/usr/bin/env python3
import argparse

from mertk.bindcraft import prepare_job

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Prepare a pinned BindCraft job for a chosen target construct"
    )
    parser.add_argument("--target-pdb", required=True)
    parser.add_argument("--target-chain", default="A")
    parser.add_argument("--output", required=True)
    parser.add_argument("--name", default="MerTK")
    parser.add_argument("--lengths", type=int, nargs=2, default=[65, 120])
    parser.add_argument("--designs", type=int, default=10)
    parser.add_argument("--hotspots", help="Optional author-numbered patch, e.g. A25,A30-33")
    args = parser.parse_args()
    output = prepare_job(
        args.target_pdb,
        args.target_chain,
        args.output,
        args.lengths,
        args.designs,
        args.hotspots,
        args.name,
    )
    print(
        f"Prepared {output / 'run.sh'}; launch from the CUDA host after BindCraft environment setup."
    )
