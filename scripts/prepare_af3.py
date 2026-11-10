#!/usr/bin/env python3
"""Export matched AF3 jobs with shared target MSAs and no structural templates."""

import argparse
import json
from pathlib import Path

from mertk.af3 import prepare_jobs

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, default=Path("configs/af3_evaluation.json"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--target-msas", type=Path, help="JSON mapping target sequence SHA256 to cached A3M path"
    )
    parser.add_argument(
        "--single-sequence",
        action="store_true",
        help="Explicit alternative protocol with no target MSA",
    )
    parser.add_argument("--samples-per-design-seed", type=int)
    parser.add_argument("--max-complexes", type=int)
    args = parser.parse_args()
    protocol = json.loads(args.protocol.read_text())
    if args.single_sequence:
        if args.target_msas:
            parser.error("Choose either cached target MSAs or --single-sequence")
        if protocol["protocol_state"] == "locked":
            parser.error(
                "Set single_sequence in a separate locked protocol instead of overriding it"
            )
        protocol["target_msa"] = "single_sequence"
    msas = None
    if args.target_msas:
        msas = {
            key: str((args.target_msas.resolve().parent / value).resolve())
            for key, value in json.loads(args.target_msas.read_text()).items()
        }
    plan = prepare_jobs(
        args.candidates,
        protocol,
        args.output,
        msas,
        args.samples_per_design_seed,
        args.max_complexes,
    )
    count = len(plan["jobs"]) * len(protocol["model_seeds"]) * protocol["samples_per_seed"]
    print(
        f"Prepared {len(plan['jobs'])} candidates / {count} AF3 predictions. No inference launched."
    )
