#!/usr/bin/env python3
"""Submit a validated, budget-bounded native-context RFD4 repair pilot."""

import argparse
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from submit_bindcraft_round import remote
from submit_protector_round import accounting

from mertk.provenance import file_hash
from mertk.target_policy import check_manifest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("manifest", type=Path)
    p.add_argument("--submit", action="store_true")
    a = p.parse_args()
    path = a.manifest.resolve()
    suite = path.parent
    r = check_manifest(path)
    if (suite / "input_validation.json").exists():
        check_manifest(path, generation=True)
    if (
        r["workflow"] != "rfd4_native_repair_v1"
        or "generate" in r["jobs"]
        or r.get("gpu_submission_intent")
    ):
        raise ValueError("Fresh validated repair submission required; reconcile previous intents")
    for name, expected in r["input_sha256"].items():
        assert file_hash(suite / name) == expected, name
    jid = r["jobs"]["author"]["job_id"]
    status = remote(["sacct", "-X", "-nP", "-j", jid, "-o", "JobID,State,ExitCode"])
    author_rows = [row.split("|") for row in status.splitlines() if row.split("|")[0] == jid]
    if len(author_rows) != 1 or author_rows[0][1] not in {"PENDING", "RUNNING", "COMPLETED"}:
        raise ValueError("Author job failed or is missing; cannot queue GPU dependency")
    if author_rows[0][1] == "COMPLETED" and author_rows[0][2] != "0:0":
        raise ValueError("Author job did not exit successfully")
    raw, rows = accounting(r)
    used = sum(row["spent_or_reserved_gpu_seconds"] for row in rows) / 3600
    budget = {
        "checked_at": datetime.now(timezone.utc).isoformat(),
        "used_or_reserved_gpu_hours": used,
        "new_gpu_hours_max": 0.5,
        "total_max_gpu_hours": used + 0.5,
        "cap": 24,
    }
    if used + 0.5 > 24:
        raise ValueError("Insufficient GPU budget")
    print(json.dumps(budget, indent=2), flush=True)
    if not a.submit:
        return
    r["budget_at_submission"] = budget
    r["gpu_submission_intent"] = {"created_at": budget["checked_at"], "reconciled": False}
    path.write_text(json.dumps(r, indent=2) + "\n")
    rem = r["remote_campaign"]
    project = r["remote_project"]
    subprocess.run(
        [
            "rsync",
            "-az",
            "-e",
            "ssh -o HostName=hoffman2.oarc.ucla.edu -o ControlMaster=no -o ControlPath=none",
            str(suite) + "/",
            "hoffman2:" + rem + "/",
        ],
        check=True,
    )
    common = [
        "sbatch",
        "--parsable",
        "--account=pi_pallard",
        "--partition=campus",
        "--qos=campus24",
        "--cpus-per-task=4",
    ]
    command = common + [
        "--mem=96G",
        "--gres=gpu:1",
        "--constraint=(A100&vmem80G)|H200",
        "--time=00:30:00",
        "--dependency=afterok:" + jid,
        "--kill-on-invalid-dep=yes",
        "--job-name=" + r["campaign_id"] + "_generate",
        "--output=" + rem + "/generate.log",
        rem + "/workflow/run_rfd4_repair.sh",
        "generate",
        project,
        rem,
    ]
    job = remote(command).strip().split(";")[0]
    if not job.isdecimal():
        raise ValueError("Unrecognized sbatch result; reconcile intent")
    r["jobs"]["generate"] = {"job_id": job, "command": command}
    r["gpu_submission_intent"]["reconciled"] = True
    path.write_text(json.dumps(r, indent=2) + "\n")
    command = common + [
        "--mem=24G",
        "--time=00:10:00",
        "--dependency=afterany:" + job,
        "--kill-on-invalid-dep=yes",
        "--job-name=" + r["campaign_id"] + "_quality",
        "--output=" + rem + "/quality.log",
        rem + "/workflow/run_rfd4_repair.sh",
        "quality",
        project,
        rem,
    ]
    job = remote(command).strip().split(";")[0]
    if not job.isdecimal():
        raise ValueError("Unrecognized QC submission result")
    r["jobs"]["quality"] = {"job_id": job, "command": command}
    path.write_text(json.dumps(r, indent=2) + "\n")
    subprocess.run(
        [
            "scp",
            "-o",
            "HostName=hoffman2.oarc.ucla.edu",
            str(path),
            "hoffman2:" + rem + "/manifest.json",
        ],
        check=True,
        capture_output=True,
    )
    print(json.dumps(r["jobs"], indent=2))


if __name__ == "__main__":
    main()
