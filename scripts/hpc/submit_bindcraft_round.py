#!/usr/bin/env python3
"""Budget-check, stage and submit a frozen BindCraft round with durable job intents."""

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
    "ControlMaster=no",
    "-o",
    "ControlPath=none",
    "-o",
    "HostName=hoffman2.oarc.ucla.edu",
    "-o",
    "StrictHostKeyChecking=yes",
]
TERMINAL = {
    "COMPLETED",
    "FAILED",
    "CANCELLED",
    "TIMEOUT",
    "OUT_OF_MEMORY",
    "NODE_FAIL",
    "PREEMPTED",
    "BOOT_FAIL",
    "DEADLINE",
    "REVOKED",
}
PROJECT_JOB_PREFIXES = ("mertk_", "course_", "round_")
ALLOCATION_ID = re.compile(r"\d+(?:[_+]\d+)?")


def seconds(value):
    h, m, s = map(int, value.split(":"))
    return 3600 * h + 60 * m + s


def gpu_count(tres):
    generic = re.search(r"(?:^|,)gres/gpu=(\d+)(?:,|$)", tres)
    if generic:
        return int(generic[1])
    return sum(int(n) for n in re.findall(r"(?:^|,)gres/gpu:[^=,]+=(\d+)", tres))


def parse_accounting(text):
    rows = []
    seen = set()
    for line in text.splitlines():
        cells = line.split("|")
        if len(cells) < 8:
            continue
        job_id, name, state, exit_code, elapsed, limit, allocated, requested = cells[:8]
        if not name.startswith(PROJECT_JOB_PREFIXES):
            continue
        if not ALLOCATION_ID.fullmatch(job_id):
            if gpu_count(allocated) or gpu_count(requested):
                raise ValueError(f"Cannot account for project GPU allocation {job_id!r}")
            continue
        if job_id in seen:
            raise ValueError("Duplicate allocation record")
        seen.add(job_id)
        gpus = gpu_count(allocated) or gpu_count(requested)
        if not gpus:
            continue
        terminal = state.split()[0] in TERMINAL
        actual = int(elapsed) * gpus
        charged = actual if terminal else max(actual, int(limit) * 60 * gpus)
        rows.append(
            {
                "job_id": job_id,
                "name": name,
                "state": state,
                "exit_code": exit_code,
                "gpus": gpus,
                "elapsed_gpu_seconds": actual,
                "spent_or_reserved_gpu_seconds": charged,
            }
        )
    return rows


def remote(command):
    return subprocess.run(
        SSH + ["hoffman2", shlex.join(command)],
        check=True,
        text=True,
        capture_output=True,
        timeout=90,
    ).stdout


def accounting(record):
    """Count each allocation once and reconcile the live queue with Slurm accounting."""
    fields = "JobID%64,JobName%100,State,ExitCode,ElapsedRaw,TimelimitRaw,AllocTRES%150,ReqTRES%150"
    query = ["sacct", "--allocations", "--array", "-nP", "-o", fields]
    user = record.get("accounting_user", "sujit009")
    raw = remote([*query, "-u", user, "-S", record["accounting_since"]])
    active = remote(["squeue", "--array", "-u", user, "-h", "-o", "%i|%j"])
    ids = []
    for line in active.splitlines():
        job_id, _, name = line.partition("|")
        if name.startswith(PROJECT_JOB_PREFIXES):
            if not ALLOCATION_ID.fullmatch(job_id):
                raise ValueError(f"Cannot reconcile live project job {job_id!r}")
            ids.append(job_id)
    if ids:
        raw += "\n" + remote([*query, "-j", ",".join(ids)])
    # The explicit job query includes held submissions that a date query may omit.
    unique = {line.split("|", 1)[0]: line for line in raw.splitlines() if "|" in line}
    if missing := set(ids) - unique.keys():
        raise ValueError(
            f"Missing accounting rows for live project jobs: {', '.join(sorted(missing))}"
        )
    rows = parse_accounting("\n".join(unique.values()))
    if not rows:
        raise ValueError("Missing project GPU accounting")
    return raw, rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--submit", action="store_true", help="Without this flag, print a dry run")
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--validation-only", action="store_true")
    modes.add_argument("--after-validation", action="store_true")
    args = parser.parse_args()
    manifest_path = args.manifest.resolve()
    record = check_manifest(manifest_path, generation=args.after_validation)
    if record.get("workflow") == "ig_followup_v1":
        raise SystemExit(
            "Use submit_ig_followup.py for the glycan-preflight and matched AF3 workflow"
        )
    if args.after_validation:
        if set(record["jobs"]) != {"input_check"} or set(record["submission_intents"]) != {
            "input_check"
        }:
            raise SystemExit(
                "Expected exactly one validation job; reconcile any previous submissions"
            )
        if not record["submission_intents"]["input_check"]["reconciled"]:
            raise SystemExit("Unreconciled input-check submission")
        job_id = record["jobs"]["input_check"]["job_id"]
        stdout = remote(
            ["sacct", "--allocations", "-nP", "-j", job_id, "-o", "JobID,State,ExitCode"]
        )
        rows = [line.split("|") for line in stdout.splitlines() if line.split("|")[0] == job_id]
        if len(rows) != 1 or rows[0][1:3] != ["COMPLETED", "0:0"]:
            raise SystemExit("Input validation must complete successfully before GPU submission")
    elif record["jobs"] or record["submission_intents"]:
        raise SystemExit("Submission already attempted; reconcile instead of blindly retrying")
    accounting_raw, accounting_rows = accounting(record)
    committed = sum(r["spent_or_reserved_gpu_seconds"] for r in accounting_rows) / 3600
    reservation = sum(seconds(s["walltime"]) for s in record["experiments"].values()) / 3600
    total = committed + reservation + record["af3_reserved_gpu_hours"]
    if total > record["course_gpu_hour_cap"]:
        raise SystemExit(
            f"Budget exceeded: {committed:.3f} committed + {reservation} new + AF3 reserve = {total:.3f} h"
        )
    budget = {
        "checked_at_utc": datetime.now(timezone.utc).isoformat(),
        "accounting": accounting_rows,
        "spent_or_reserved_gpu_hours": committed,
        "new_maximum_gpu_hours": reservation,
        "af3_reserved_gpu_hours": record["af3_reserved_gpu_hours"],
        "total_with_af3_reserve": total,
        "cap": record["course_gpu_hour_cap"],
    }
    print(json.dumps({k: v for k, v in budget.items() if k != "accounting"}, indent=2))
    specs = [("input_check", "validate", "", None, False, "00:10:00")]
    if args.after_validation:
        specs = [
            (
                "pdl1_control_generate",
                "generate",
                "pdl1_control",
                "input_check",
                True,
                record["experiments"]["pdl1_control"]["walltime"],
            ),
            ("control_gate", "control-gate", "", "any:pdl1_control_generate", False, "00:10:00"),
        ]
        for key in record["experiments"]:
            if key == "pdl1_control":
                continue
            specs.extend(
                [
                    (
                        f"{key}_generate",
                        "generate",
                        key,
                        "control_gate",
                        True,
                        record["experiments"][key]["walltime"],
                    ),
                    (f"{key}_quality", "report", key, f"any:{key}_generate", False, "00:20:00"),
                ]
            )
    if not args.submit:
        print(json.dumps(specs, indent=2))
        return
    record["budget_at_submission"] = budget
    manifest_path.write_text(json.dumps(record, indent=2) + "\n")
    remote(["mkdir", "-p", record["remote_campaign"]])
    subprocess.run(
        [
            "rsync",
            "-az",
            "--exclude=__pycache__",
            "-e",
            shlex.join(SSH),
            str(manifest_path.parent) + "/",
            f"hoffman2:{record['remote_campaign']}/",
        ],
        check=True,
        timeout=180,
    )
    for name, stage, experiment, dependency, gpu, limit in specs:
        command = [
            "sbatch",
            "--parsable",
            f"--account={record['authorized_account']}",
            f"--partition={record['partition']}",
            f"--qos={record['qos']}",
            "--nodes=1",
            "--cpus-per-task=8",
            f"--mem={'64G' if gpu else '16G'}",
            f"--time={limit}",
            f"--job-name={record['campaign_id']}_{name}",
            f"--output={record['remote_project']}/logs/%x_%j.log",
            f"--chdir={record['remote_project']}",
        ]
        if gpu:
            command += ["--gres=gpu:1", f"--constraint={record.get('gpu_constraint', 'A100')}"]
        if dependency:
            any_state = dependency.startswith("any:")
            parent = dependency.split(":")[-1]
            command += [
                f"--dependency={'afterany' if any_state else 'afterok'}:{record['jobs'][parent]['job_id']}",
                "--kill-on-invalid-dep=yes",
            ]
        command += [
            f"{record['remote_campaign']}/workflow/run_bindcraft_round.sh",
            stage,
            record["remote_project"],
            record["remote_campaign"],
        ]
        if experiment:
            command.append(experiment)
        record["submission_intents"][name] = {
            "command": command,
            "reconciled": False,
            "requested_at_utc": datetime.now(timezone.utc).isoformat(),
        }
        manifest_path.write_text(json.dumps(record, indent=2) + "\n")
        stdout = remote(command)
        ids = [
            line.split(";")[0]
            for line in stdout.splitlines()
            if re.fullmatch(r"\d+(;[\w-]+)?", line)
        ]
        if len(ids) != 1:
            raise SystemExit(
                "Uncertain submission response; reconcile the recorded job name before retrying"
            )
        record["jobs"][name] = {"job_id": ids[0], "command": command}
        record["submission_intents"][name]["reconciled"] = True
        manifest_path.write_text(json.dumps(record, indent=2) + "\n")
        print(f"{name}: {ids[0]}", flush=True)
    subprocess.run(
        [
            "rsync",
            "-az",
            "-e",
            shlex.join(SSH),
            str(manifest_path),
            f"hoffman2:{record['remote_campaign']}/manifest.json",
        ],
        check=True,
        timeout=90,
    )


if __name__ == "__main__":
    main()
