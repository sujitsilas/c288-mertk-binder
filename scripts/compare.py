#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

from mertk.comparison import run_comparison

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Evaluate a matched sequence-design comparison arm"
    )
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--arm", required=True)
    parser.add_argument("--checkpoint", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--samples", type=int, default=8)
    parser.add_argument("--sampling-steps", type=int, default=16)
    parser.add_argument("--temperature", type=float, default=1.0)
    parser.add_argument("--split", choices=["val", "test", "application_holdout"], default="test")
    parser.add_argument("--expected-dataset")
    parser.add_argument("--application-manifest", type=Path)
    parser.add_argument("--smoke-test", action="store_true")
    args = parser.parse_args()
    result = run_comparison(
        json.loads(args.config.read_text()),
        args.arm,
        args.output,
        args.checkpoint,
        args.samples,
        args.sampling_steps,
        args.temperature,
        args.split,
        args.smoke_test,
        args.expected_dataset,
        args.application_manifest,
    )
    print(f"{result['arm']}: {len(result['cases'])} matched examples")
