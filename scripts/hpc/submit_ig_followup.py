#!/usr/bin/env python3
"""Submit the bounded Ig/glycan/AF3 DAG with durable, reconciled job intents."""

import argparse
import json
import re
import shlex
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from bindcraft_round import validate
from submit_bindcraft_round import SSH, accounting, remote, seconds

from mertk.target_policy import check_manifest


def job_specs(record):
    specs = [
        ("input_check", "validate", "", [], "afterok", False, "00:10:00"),
        ("glycan_preflight", "glycan-preflight", "", ["input_check"], "afterok", True, "00:20:00"),
    ]
    for key, spec in record["experiments"].items():
        specs.extend(
            [
                (
                    f"{key}_generate",
                    "generate",
                    key,
                    ["glycan_preflight"],
                    "afterok",
                    True,
                    spec["walltime"],
                ),
                (
                    f"{key}_quality",
                    "quality",
                    key,
                    [f"{key}_generate"],
                    "afterany",
                    False,
                    "00:15:00",
                ),
            ]
        )
    quality = [f"{key}_quality" for key in record["experiments"]]
    specs.extend(
        [
            ("protector_triage", "protector-triage", "", quality, "afterok", False, "00:10:00"),
            ("af3_parents", "af3", "parents", ["protector_triage"], "afterok", True, "00:45:00"),
            ("af3_tulp2", "af3", "mertk_tulp2", ["input_check"], "afterok", True, "00:45:00"),
            ("af3_mcherry", "af3", "mertk_mcherry", ["input_check"], "afterok", True, "00:30:00"),
            (
                "af3_report",
                "report",
                "",
                ["af3_parents", "af3_tulp2", "af3_mcherry"],
                "afterany",
                False,
                "00:15:00",
            ),
        ]
    )
    return specs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--submit", action="store_true")
    args = parser.parse_args()
    path = args.manifest.resolve()
    record = check_manifest(path, generation=True)
    validate(path.parent)
    if record.get("workflow") != "ig_followup_v1" or record["jobs"] or record["submission_intents"]:
        raise ValueError("Requires a new prepared Ig follow-up; reconcile any existing submission")
    raw, accounting_rows = accounting(record)
    specs = job_specs(record)
    committed = sum(r["spent_or_reserved_gpu_seconds"] for r in accounting_rows) / 3600
    new = sum(seconds(limit) for _, _, _, _, _, gpu, limit in specs if gpu) / 3600
    if committed + new > record["course_gpu_hour_cap"]:
        raise ValueError(f"Budget exceeded: {committed:.3f} + {new:.3f} hours")
    budget = {
        "checked_at": datetime.now(timezone.utc).isoformat(),
        "spent_or_reserved_gpu_hours": committed,
        "new_maximum_gpu_hours": new,
        "total_gpu_hours": committed + new,
        "cap": record["course_gpu_hour_cap"],
        "af3_included_hours": 2.0,
        "accounting": accounting_rows,
    }
    print(json.dumps({k: v for k, v in budget.items() if k != "accounting"}, indent=2), flush=True)
    if not args.submit:
        print(json.dumps(specs, indent=2))
        return
    record["budget_at_submission"] = budget
    path.write_text(json.dumps(record, indent=2) + "\n")
    (path.parent / "accounting_at_submission.txt").write_text(raw)
    remote(["mkdir", "-p", record["remote_campaign"]])
    subprocess.run(
        [
            "rsync",
            "-az",
            "--exclude=__pycache__",
            "-e",
            shlex.join(SSH),
            str(path.parent) + "/",
            f"hoffman2:{record['remote_campaign']}/",
        ],
        check=True,
        timeout=180,
    )
    for name, stage, experiment, parents, kind, gpu, limit in specs:
        command = [
            "sbatch",
            "--parsable",
            f"--account={record['authorized_account']}",
            f"--partition={record['partition']}",
            f"--qos={record['qos']}",
            "--nodes=1",
            f"--cpus-per-task={record.get('af3_cpus_per_task', 8) if stage == 'af3' else 8}",
            f"--mem={'96G' if stage == 'af3' else '64G' if gpu else '16G'}",
            f"--time={limit}",
            f"--job-name={record['campaign_id']}_{name}",
            f"--output={record['remote_project']}/logs/%x_%j.log",
            f"--chdir={record['remote_project']}",
        ]
        if gpu:
            constraint = (
                record["af3_gpu_constraint"] if stage == "af3" else record["gpu_constraint"]
            )
            command += ["--gres=gpu:1", f"--constraint={constraint}"]
        if parents:
            ids = ":".join(record["jobs"][p]["job_id"] for p in parents)
            command += [f"--dependency={kind}:{ids}", "--kill-on-invalid-dep=yes"]
        command += [
            f"{record['remote_campaign']}/workflow/run_ig_followup.sh",
            stage,
            record["remote_project"],
            record["remote_campaign"],
        ]
        if experiment:
            command.append(experiment)
        record["submission_intents"][name] = {
            "command": command,
            "reconciled": False,
            "requested_at": datetime.now(timezone.utc).isoformat(),
        }
        path.write_text(json.dumps(record, indent=2) + "\n")
        response = remote(command)
        ids = [
            line.split(";")[0]
            for line in response.splitlines()
            if re.fullmatch(r"\d+(;[\w-]+)?", line)
        ]
        if len(ids) != 1:
            raise RuntimeError("Uncertain submission; reconcile the named intent before retrying")
        record["jobs"][name] = {"job_id": ids[0], "command": command}
        record["submission_intents"][name]["reconciled"] = True
        path.write_text(json.dumps(record, indent=2) + "\n")
        print(f"{name}: {ids[0]}", flush=True)
    subprocess.run(
        [
            "rsync",
            "-az",
            "-e",
            shlex.join(SSH),
            str(path),
            f"hoffman2:{record['remote_campaign']}/manifest.json",
        ],
        check=True,
        timeout=90,
    )


if __name__ == "__main__":
    main()
