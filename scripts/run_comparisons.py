#!/usr/bin/env python3
"""Plan or execute matched baseline/ablation runs; subprocesses release model memory."""

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
    parser.add_argument("--base-config", type=Path, default=Path("configs/train.json"))
    parser.add_argument("--suite", type=Path, default=Path("configs/comparison_core.json"))
    parser.add_argument("--output", type=Path, default=Path("outputs/comparisons_core"))
    parser.add_argument(
        "--execute", action="store_true", help="Run jobs; default only writes a reviewable plan"
    )
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--smoke-test", action="store_true")
    parser.add_argument("--train-family-fraction", type=float)
    parser.add_argument("--data-seed", type=int)
    args = parser.parse_args()
    base, suite = json.loads(args.base_config.read_text()), json.loads(args.suite.read_text())
    for key in ("train_family_fraction", "data_seed"):
        if getattr(args, key) is not None:
            base[key] = getattr(args, key)
    if args.smoke_test:
        suite.update(seeds=suite["seeds"][:1], samples=2, sampling_steps=2)
        base.update(
            steps=2,
            accumulation_steps=1,
            eval_every=1,
            device="cpu",
            precision="fp32",
            cpu_threads=2,
        )
    identity = dataset_identity(base["manifest"], args.smoke_test)
    if not args.smoke_test:
        rows = read_records(base["manifest"])
        audit_splits(rows)
        missing = {"train", "val", suite["split"]} - {r["split"] for r in rows}
        if missing and args.execute:
            raise ValueError(
                f"Missing independent dataset splits: {sorted(missing)}; expand data first"
            )
        if missing:
            print(f"Plan only: missing dataset splits {sorted(missing)}")
    definition = {
        "runner_version": 2,
        "base": base,
        "suite": suite,
        "dataset_identity": identity,
        "smoke_test": args.smoke_test,
    }
    plan_id = hashlib.sha256(json.dumps(definition, sort_keys=True).encode()).hexdigest()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=True)
    plan_file = root / "plan.json"
    if plan_file.exists() and json.loads(plan_file.read_text())["plan_id"] != plan_id:
        raise ValueError("Existing plan differs; choose another output directory")
    jobs = []
    sampling_protocol = {
        "samples": suite["samples"],
        "steps": suite["sampling_steps"],
        "temperatures": {
            arm: suite.get("temperatures", {}).get(arm, suite["temperature"])
            for arm in suite["arms"]
        },
    }
    for arm, overrides in suite["arms"].items():
        for seed in suite["seeds"]:
            directory = root / f"{arm}_seed{seed}"
            config = dict(
                base,
                **overrides,
                seed=seed,
                output=str(directory),
                sampling_protocol=sampling_protocol,
            )
            config_file = root / f"{arm}_seed{seed}.json"
            config_file.write_text(json.dumps(config, indent=2) + "\n")
            checkpoint = directory / "best.pt"
            resume_checkpoint = directory / "last.pt"
            trainable = arm not in {"frozen_esm3", "proteinmpnn"}
            train_command = (
                [sys.executable, "scripts/train.py", "--config", str(config_file)]
                if trainable
                else None
            )
            if train_command and args.smoke_test:
                train_command += ["--smoke-test", "--output", str(directory)]
            result = directory / "comparison.json"
            compare = [
                sys.executable,
                "scripts/compare.py",
                "--config",
                str(config_file),
                "--arm",
                arm,
                "--output",
                str(result),
                "--samples",
                str(suite["samples"]),
                "--sampling-steps",
                str(suite["sampling_steps"]),
                "--temperature",
                str(sampling_protocol["temperatures"][arm]),
                "--split",
                suite["split"],
                "--expected-dataset",
                identity,
            ]
            if trainable:
                compare += ["--checkpoint", str(checkpoint)]
            if args.smoke_test:
                compare += ["--smoke-test"]
            jobs.append(
                {
                    "arm": arm,
                    "seed": seed,
                    "train": train_command,
                    "compare": compare,
                    "checkpoint": str(checkpoint),
                    "resume_checkpoint": str(resume_checkpoint),
                    "result": str(result),
                }
            )
    plan_file.write_text(
        json.dumps({"plan_id": plan_id, **definition, "jobs": jobs}, indent=2) + "\n"
    )
    print(f"Prepared {len(jobs)} arm/seed jobs in {plan_file}", flush=True)
    if not args.execute:
        return
    for job in jobs:
        if Path(job["result"]).exists():
            if not args.resume:
                raise FileExistsError(f"{job['result']} already exists; use --resume")
            continue
        print(f"Running {job['arm']} seed {job['seed']}", flush=True)
        if job["train"]:
            command = job["train"]
            if args.resume and Path(job["resume_checkpoint"]).exists():
                command += ["--resume", job["resume_checkpoint"]]
            subprocess.run(command, check=True)
        subprocess.run(job["compare"], check=True)
    results = [json.loads(Path(job["result"]).read_text()) for job in jobs]
    report = summarize(results)
    (root / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    if not args.smoke_test:
        export_candidates(read_records(base["manifest"]), results, root / "candidates.jsonl")
    print(f"Matched family-level summary: {root / 'summary.json'}")


if __name__ == "__main__":
    main()
