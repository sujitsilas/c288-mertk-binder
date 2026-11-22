#!/usr/bin/env python3
"""Synchronize a bounded RFD4 repair campaign without submitting new work."""

import argparse
import json
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

from submit_bindcraft_round import TERMINAL, remote


def refresh(path):
    suite = path.parent
    r = json.loads(path.read_text())
    ids = ",".join(row["job_id"] for row in r["jobs"].values())
    raw = remote(["sacct", "-X", "-nP", "-j", ids, "-o", "JobID,State,ExitCode,Elapsed,NodeList"])
    rows = {
        line.split("|")[0]: line.split("|")[1:5]
        for line in raw.splitlines()
        if line.split("|")[0].isdigit()
    }
    subprocess.run(
        [
            "rsync",
            "-az",
            "--exclude=manifest.json",
            "-e",
            "ssh -o HostName=hoffman2.oarc.ucla.edu -o ControlMaster=no -o ControlPath=none",
            "hoffman2:" + r["remote_campaign"] + "/",
            str(suite) + "/",
        ],
        check=True,
        capture_output=True,
    )
    jobs = {
        key: dict(
            job_id=value["job_id"],
            **dict(
                zip(
                    ["state", "exit_code", "elapsed", "nodes"],
                    rows.get(value["job_id"], ["UNKNOWN", "", "", ""]),
                )
            ),
        )
        for key, value in r["jobs"].items()
    }
    done = all(v["state"].split()[0] in TERMINAL for v in jobs.values())
    result = {
        "checked_at": datetime.now(timezone.utc).isoformat(),
        "jobs": jobs,
        "all_jobs_terminal": done,
    }
    (suite / "status.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result), flush=True)
    subprocess.run(
        [
            str(Path(".venv/bin/python").absolute()),
            "scripts/build_prediction_review.py",
            "outputs/hpc/round_20261114_ig_guard_01",
        ],
        check=True,
        capture_output=True,
    )
    return done


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("manifest", type=Path)
    p.add_argument("--watch", action="store_true")
    a = p.parse_args()
    while True:
        try:
            done = refresh(a.manifest.resolve())
            if done or not a.watch:
                break
        except Exception as error:
            print(str(error), flush=True)
            if not a.watch:
                raise
        time.sleep(120)
