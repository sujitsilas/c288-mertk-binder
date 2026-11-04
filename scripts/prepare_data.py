#!/usr/bin/env python3
import argparse

from mertk.data import build_dataset

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract leakage-resistant chain-pair examples")
    parser.add_argument("--manifest", default="data/raw/manifest.json")
    parser.add_argument("--output", default="data/processed")
    parser.add_argument("--max-length", type=int, default=512)
    parser.add_argument("--min-interface", type=int, default=5)
    args = parser.parse_args()
    records = build_dataset(args.manifest, args.output, args.max_length, args.min_interface)
    print(f"Prepared {len(records)} directed interfaces; cluster chains before training.")
    if not records:
        raise SystemExit("No usable interfaces. Expand the download or inspect the length filter.")
