#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

from mertk.structural import score_bundle

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Score every planned candidate with matched AF3 structural gates"
    )
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(score_bundle(args.manifest, args.predictions, args.output), indent=2))
