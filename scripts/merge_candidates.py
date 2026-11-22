#!/usr/bin/env python3
"""Combine candidate manifests after validating target and epitope coverage."""

import argparse
import json
from pathlib import Path

from mertk.af3 import checked_reference, validate_cohort

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", nargs="+", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, default=Path("configs/af3_evaluation.json"))
    args = parser.parse_args()
    rows = [
        json.loads(line)
        for path in args.inputs
        for line in path.read_text().splitlines()
        if line.strip()
    ]
    validate_cohort(rows, json.loads(args.protocol.read_text())["required_arms"])
    for row in rows:
        checked_reference(row)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as handle:
        handle.write("".join(json.dumps(row) + "\n" for row in rows))
    print(f"Combined {len(rows)} candidates")
