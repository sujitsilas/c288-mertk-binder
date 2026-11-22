#!/usr/bin/env python3
"""Submit an independently budgeted stalk-protector DAG; hold AF3 until eligible input exists."""

import argparse
import json
import re
import shlex
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from bindcraft_round import validate
from submit_bindcraft_round import SSH, remote, seconds
from submit_bindcraft_round import accounting as accounting

from mertk.target_policy import check_manifest


def job_specs(record):
    rows = [
        ("input_check", "validate", "", [], "afterok", False, "00:10:00", False),
        (
            "protector_preflight",
            "preflight",
            "",
            ["input_check"],
            "afterok",
            True,
            "00:10:00",
            False,
        ),
    ]
    for key, spec in record["experiments"].items():
        rows += [
            (
                key + "_generate",
                "generate",
                key,
                ["protector_preflight"],
                "afterok",
                True,
                spec["walltime"],
                False,
            ),
            (
                key + "_quality",
                "quality",
                key,
                [key + "_generate"],
                "afterany",
                False,
                "00:15:00",
                False,
            ),
        ]
    rows += [
        (
            "protector_select",
            "select",
            "",
            [k + "_quality" for k in record["experiments"]],
            "afterany",
            False,
            "00:10:00",
            False,
        ),
        ("af3_assessment", "af3", "", ["protector_select"], "afterok", True, "00:20:00", True),
        (
            "protector_report",
            "report",
            "",
            ["protector_select", "af3_assessment"],
            "afterany",
            False,
            "00:15:00",
            False,
        ),
    ]
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--submit", action="store_true")
    args = parser.parse_args()
    path = args.manifest.resolve()
    record = check_manifest(path, generation=True)
    validate(path.parent)
    if (
        record.get("workflow") != "cleavage_protector_v1"
        or record["jobs"]
        or record["submission_intents"]
    ):
        raise ValueError("Requires fresh dedicated protector suite; reconcile prior submissions")
    raw, rows = accounting(record)
    specs = job_specs(record)
    committed = sum(r["spent_or_reserved_gpu_seconds"] for r in rows) / 3600
    reserved = sum(seconds(r[6]) for r in specs if r[5]) / 3600
    if committed + reserved > record["course_gpu_hour_cap"]:
        raise ValueError(f"GPU budget exceeded: {committed} + {reserved}")
    budget = dict(
        checked_at=datetime.now(timezone.utc).isoformat(),
        spent_or_reserved_gpu_hours=committed,
        new_maximum_gpu_hours=reserved,
        total_gpu_hours=committed + reserved,
        cap=record["course_gpu_hour_cap"],
        accounting=rows,
    )
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
    for name, stage, experiment, parents, dependency, gpu, limit, held in specs:
        cmd = [
            "sbatch",
            "--parsable",
            f"--account={record['authorized_account']}",
            "--partition=campus",
            "--qos=campus24",
            "--nodes=1",
            "--cpus-per-task=4",
            f"--mem={'96G' if stage == 'af3' else '64G' if gpu else '16G'}",
            f"--time={limit}",
            f"--job-name={record['campaign_id']}_{name}",
            f"--output={record['remote_project']}/logs/%x_%j.log",
            f"--chdir={record['remote_project']}",
        ]
        if gpu:
            cmd += [
                "--gres=gpu:1",
                "--constraint="
                + record["af3_gpu_constraint" if stage == "af3" else "gpu_constraint"],
            ]
        if parents:
            cmd += [
                f"--dependency={dependency}:"
                + ":".join(record["jobs"][p]["job_id"] for p in parents),
                "--kill-on-invalid-dep=yes",
            ]
        if held:
            cmd += ["--hold"]
        cmd += [
            record["remote_campaign"] + "/workflow/run_protector_round.sh",
            stage,
            record["remote_project"],
            record["remote_campaign"],
        ]
        if experiment:
            cmd.append(experiment)
        record["submission_intents"][name] = dict(
            command=cmd, reconciled=False, requested_at=datetime.now(timezone.utc).isoformat()
        )
        path.write_text(json.dumps(record, indent=2) + "\n")
        reply = remote(cmd)
        ids = [x.split(";")[0] for x in reply.splitlines() if re.fullmatch(r"\d+(;[\w-]+)?", x)]
        if len(ids) != 1:
            raise RuntimeError("Uncertain submission; reconcile this job name before retrying")
        record["jobs"][name] = dict(job_id=ids[0], command=cmd)
        record["submission_intents"][name]["reconciled"] = True
        path.write_text(json.dumps(record, indent=2) + "\n")
        print(name, ids[0], flush=True)
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
