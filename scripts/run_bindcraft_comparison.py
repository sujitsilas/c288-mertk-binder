#!/usr/bin/env python3
"""Redesign imported BindCraft holdouts with an existing trained comparison suite."""

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

from mertk.candidates import export_candidates
from mertk.comparison import dataset_identity, summarize
from mertk.dataset import read_records
from mertk.training import audit_splits


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--trained-suite", type=Path, required=True, help="Training suite plan.json"
    )
    parser.add_argument("--manifest", type=Path, default=Path("data/application/examples.jsonl"))
    parser.add_argument("--output", type=Path, default=Path("outputs/bindcraft_comparison"))
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    trained = json.loads(args.trained_suite.read_text())
    if trained["smoke_test"]:
        raise ValueError("Train on real family-split data before application redesign")
    base, suite = trained["base"], trained["suite"]
    if dataset_identity(base["manifest"]) != trained["dataset_identity"]:
        raise ValueError("Training dataset differs from the saved suite plan")
    rows = read_records(args.manifest)
    if any(row["split"] != "application_holdout" or "source_pdb" not in row for row in rows):
        raise ValueError("Import accepted complexes with scripts/import_bindcraft.py first")
    audit_splits(read_records(base["manifest"]) + rows)
    identity = dataset_identity(args.manifest)
    root = args.output.resolve()
    definition = {
        "runner_version": 2,
        "training_plan_id": trained["plan_id"],
        "dataset_identity": identity,
        "manifest": str(args.manifest.resolve()),
    }
    plan_id = hashlib.sha256(json.dumps(definition, sort_keys=True).encode()).hexdigest()
    root.mkdir(parents=True, exist_ok=True)
    plan_file = root / "plan.json"
    if plan_file.exists() and json.loads(plan_file.read_text())["plan_id"] != plan_id:
        raise ValueError("Existing application plan differs; choose another output directory")
    jobs = []
    sampling_protocol = {
        "samples": suite["samples"],
        "steps": suite["sampling_steps"],
        "temperatures": {
            arm: suite.get("temperatures", {}).get(arm, suite["temperature"])
            for arm in suite["arms"]
        },
    }
    for trained_job in trained["jobs"]:
        arm, seed = trained_job["arm"], trained_job["seed"]
        config = dict(base, **suite["arms"][arm], seed=seed, sampling_protocol=sampling_protocol)
        config_file = root / f"{arm}_seed{seed}.json"
        config_file.write_text(json.dumps(config, indent=2) + "\n")
        result = root / f"{arm}_seed{seed}" / "comparison.json"
        command = [
            sys.executable,
            "scripts/compare.py",
            "--config",
            str(config_file),
            "--arm",
            arm,
            "--output",
            str(result),
            "--split",
            "application_holdout",
            "--application-manifest",
            str(args.manifest.resolve()),
            "--expected-dataset",
            identity,
            "--samples",
            str(suite["samples"]),
            "--sampling-steps",
            str(suite["sampling_steps"]),
            "--temperature",
            str(sampling_protocol["temperatures"][arm]),
        ]
        if trained_job["train"]:
            checkpoint = Path(trained_job["checkpoint"])
            if args.execute and not checkpoint.is_file():
                raise FileNotFoundError(f"Complete the training comparison first: {checkpoint}")
            command += ["--checkpoint", str(checkpoint)]
        jobs.append({"arm": arm, "seed": seed, "command": command, "result": str(result)})
    plan_file.write_text(
        json.dumps({"plan_id": plan_id, **definition, "jobs": jobs}, indent=2) + "\n"
    )
    print(f"Prepared {len(jobs)} redesign jobs; original BindCraft sequences retained as controls.")
    if not args.execute:
        return
    for job in jobs:
        if Path(job["result"]).exists():
            if not args.resume:
                raise FileExistsError(f"{job['result']} exists; use --resume")
            continue
        subprocess.run(job["command"], check=True)
    results = [json.loads(Path(job["result"]).read_text()) for job in jobs]
    report = summarize(results)
    (root / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    export_candidates(rows, results, root / "candidates.jsonl")
    print(
        f"Saved seed-retention/diversity summary and candidates for structure evaluation in {root}"
    )


if __name__ == "__main__":
    main()
