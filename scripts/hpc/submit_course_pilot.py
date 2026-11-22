#!/usr/bin/env python3
"""Submit bounded jobs against existing runtimes; retain intents to avoid blind retries."""

import argparse
import json
import re
import shlex
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from mertk.target_policy import check_manifest

SSH = [
    "ssh",
    "-o",
    "BatchMode=yes",
    "-o",
    "ConnectTimeout=15",
    "-o",
    "HostName=hoffman2.oarc.ucla.edu",
    "-o",
    "UserKnownHostsFile=outputs/hpc/known_hosts",
    "hoffman2",
]


def validation_succeeded(stdout, job_id):
    """Accept parsable sacct records with or without their optional final delimiter."""
    rows = [
        line.strip().split("|")
        for line in stdout.splitlines()
        if line.strip().split("|")[0] == job_id
    ]
    return len(rows) == 1 and rows[0][1:3] == ["COMPLETED", "0:0"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--validation-only", action="store_true")
    mode.add_argument("--after-validation", action="store_true")
    args = parser.parse_args()
    record = check_manifest(args.manifest, generation=not args.validation_only)
    if args.after_validation:
        if (
            set(record.get("jobs", {})) != {"input_check"}
            or set(record.get("submission_intents", {})) != {"input_check"}
            or not record["submission_intents"]["input_check"]["reconciled"]
        ):
            raise SystemExit("Expected exactly one reconciled validation submission")
        job_id = record["jobs"]["input_check"]["job_id"]
        if not job_id.isdecimal():
            raise ValueError("Invalid validation job ID")
        result = subprocess.run(
            SSH
            + [
                shlex.join(
                    ["sacct", "--allocations", "-nP", "-j", job_id, "-o", "JobID,State,ExitCode"]
                )
            ],
            check=True,
            text=True,
            capture_output=True,
            timeout=90,
        )
        if not validation_succeeded(result.stdout, job_id):
            raise SystemExit("Validation is not recorded as COMPLETED with exit 0:0")
    elif record.get("jobs") or record.get("submission_intents"):
        raise SystemExit(
            "Submission already attempted. Reconcile existing jobs; do not submit duplicates."
        )
    project, campaign = record["remote_project"], record["remote_campaign"]
    specs = [
        ("input_check", "validate", None, False, "00:20:00", "16G"),
        ("bindcraft_generate", "bindcraft", "input_check", True, "04:00:00", "64G"),
        ("rfd4_control", "rfd4_control", "input_check", True, "01:00:00", "64G"),
        ("rfd4_proteina_generate", "rfd4_proteina", "rfd4_control", True, "02:00:00", "64G"),
        ("native_quality", "quality", "afterany", False, "00:20:00", "16G"),
    ]
    if args.validation_only:
        specs = specs[:1]
    elif args.after_validation:
        specs = specs[1:]
    record.setdefault("submission_intents", {})
    for name, stage, dependency, gpu, limit, memory in specs:
        job_name = f"{record['campaign_id']}_{name}"
        command = [
            "sbatch",
            "--parsable",
            f"--account={record['authorized_account']}",
            f"--partition={record['partition']}",
            f"--qos={record['qos']}",
            "--nodes=1",
            "--cpus-per-task=8",
            f"--mem={memory}",
            f"--time={limit}",
            f"--job-name={job_name}",
            f"--output={project}/logs/%x_%j.log",
            f"--chdir={project}",
        ]
        if gpu:
            command += ["--gres=gpu:1", "--constraint=A100"]
        if dependency == "afterany":
            ids = [
                record["jobs"][key]["job_id"]
                for key in ["bindcraft_generate", "rfd4_proteina_generate"]
            ]
            command += ["--dependency=afterany:" + ":".join(ids)]
        elif dependency:
            command += [
                f"--dependency=afterok:{record['jobs'][dependency]['job_id']}",
                "--kill-on-invalid-dep=yes",
            ]
        command += [f"{project}/scripts/hpc/run_course_stage.sh", stage, project, campaign]
        record["submission_intents"][name] = {
            "job_name": job_name,
            "command": command,
            "requested_at_utc": datetime.now(timezone.utc).isoformat(),
            "reconciled": False,
        }
        args.manifest.write_text(json.dumps(record, indent=2) + "\n")
        result = subprocess.run(
            SSH + [shlex.join(command)], check=True, text=True, stdout=subprocess.PIPE, timeout=90
        )
        ids = [
            line.split(";")[0]
            for line in result.stdout.splitlines()
            if re.fullmatch(r"\d+(;[A-Za-z0-9_-]+)?", line.strip())
        ]
        if len(ids) != 1:
            raise SystemExit(
                "Uncertain submission response: reconcile the recorded job name before retrying"
            )
        record["jobs"][name] = {
            "job_id": ids[0],
            "command": command,
            "submitted_at_utc": datetime.now(timezone.utc).isoformat(),
        }
        record["submission_intents"][name]["reconciled"] = True
        args.manifest.write_text(json.dumps(record, indent=2) + "\n")
        print(f"{name}: {ids[0]}", flush=True)


if __name__ == "__main__":
    main()
