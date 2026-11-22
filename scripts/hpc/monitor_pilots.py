#!/usr/bin/env python3
"""Poll Slurm, copy this campaign's logs/results, and render completed structures."""

import argparse
import json
import shlex
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

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
ROOT = Path(__file__).resolve().parents[2]
SSH = [
    "ssh",
    "-o",
    "BatchMode=yes",
    "-o",
    "ConnectTimeout=15",
    "-o",
    "ServerAliveInterval=15",
    "-o",
    "ServerAliveCountMax=2",
    "-o",
    "HostName=hoffman2.oarc.ucla.edu",
    "-o",
    f'UserKnownHostsFile="{ROOT / "outputs/hpc/known_hosts"}"',
]


def af3_status(jobs):
    """Describe observed scheduler state, not the original submission intent."""
    followup = {
        name: job.get("state", "UNKNOWN")
        for name, job in jobs.items()
        if name.startswith("af3_") and name != "af3_report"
    }
    if followup:
        return (
            "; ".join(f"{name}: {state}" for name, state in followup.items())
            + ". Inspect outputs; completion is not binding evidence."
        )
    stage = jobs.get("alphafold3_evaluate", {})
    state = stage.get("state", "NOT_SUBMITTED").split()[0]
    preparation = jobs.get("prepare_assessment", {}).get("state", "UNKNOWN").split()[0]
    if state == "COMPLETED":
        return "AF3 job completed; inspect scoring reports for prediction validity."
    if state in TERMINAL:
        return f"AF3 job {state}; candidate preparation {preparation}. No completed assessment inferred."
    if state == "RUNNING":
        return "AF3 evaluation is running."
    return f"AF3 job {state}; candidate preparation {preparation}."


def refresh_prediction_review(manifest_path, manifest):
    """Refresh the local HTML after synchronization without releasing GPU jobs."""
    config = ROOT / "configs/prediction_review.json"
    if not manifest.get("experiments") or not config.exists():
        return
    if manifest.get("workflow") == "cleavage_protector_v1":
        subprocess.run(
            [
                sys.executable,
                str(manifest_path.parent / "workflow/protector_round.py"),
                "report",
                str(manifest_path.parent),
            ],
            check=True,
            cwd=ROOT,
            capture_output=True,
            text=True,
            timeout=120,
        )
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "scripts/build_protector_review.py"),
                str(manifest_path.parent),
            ],
            check=True,
            cwd=ROOT,
            capture_output=True,
            text=True,
            timeout=120,
        )
        # Preserve the Ig viewer while refreshing its linked campaign queue.
        result = subprocess.run(
            [
                sys.executable,
                str(ROOT / "scripts/build_prediction_review.py"),
                str(ROOT / "outputs/hpc/round_20261114_ig_guard_01"),
            ],
            check=True,
            cwd=ROOT,
            capture_output=True,
            text=True,
            timeout=120,
        )
        print("Protector review refreshed; " + result.stdout.strip(), flush=True)
        return
    if manifest.get("workflow") == "ig_followup_v1":
        subprocess.run(
            [
                sys.executable,
                str(manifest_path.parent / "workflow/ig_followup.py"),
                "report",
                str(manifest_path.parent),
            ],
            check=True,
            cwd=ROOT,
            capture_output=True,
            text=True,
            timeout=120,
        )
    result = subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/build_prediction_review.py"),
            str(manifest_path.parent),
            "--config",
            str(config),
        ],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=120,
    )
    if result.returncode:
        print(
            f"Review refresh failed; scheduler monitoring continues: {result.stderr[-1500:]}",
            file=sys.stderr,
            flush=True,
        )
    else:
        print(result.stdout.strip(), flush=True)


def check(manifest_path):
    manifest = json.loads(manifest_path.read_text())
    jobs = manifest["jobs"]
    if not jobs:
        raise ValueError("No submitted jobs recorded")
    ids = [job["job_id"] for job in jobs.values()]
    if not all(job_id.isdecimal() for job_id in ids):
        raise ValueError("Invalid job identifier")
    command = [
        "sacct",
        "--allocations",
        "-nP",
        "-j",
        ",".join(ids),
        "-o",
        "JobID,State,ExitCode,Elapsed,NodeList",
    ]
    result = subprocess.run(
        SSH + ["hoffman2", shlex.join(command)],
        check=True,
        text=True,
        capture_output=True,
        timeout=90,
    )
    states = {}
    for line in result.stdout.splitlines():
        values = line.split("|")
        if len(values) >= 5 and values[0] in ids:
            states[values[0]] = dict(
                zip(("job_id", "state", "exit_code", "elapsed", "nodes"), values[:5])
            )
    snapshot = {
        "checked_at": datetime.now(timezone.utc).isoformat(),
        "jobs": {
            name: states.get(job["job_id"], {"state": "UNKNOWN"}) for name, job in jobs.items()
        },
    }
    snapshot["af3_inference"] = af3_status(snapshot["jobs"])
    done = all(
        states.get(job_id, {}).get("state", "UNKNOWN").split()[0] in TERMINAL for job_id in ids
    )
    snapshot["all_jobs_terminal"] = done
    (manifest_path.parent / "status.json").write_text(json.dumps(snapshot, indent=2) + "\n")
    for source, destination in (
        (manifest["remote_campaign"], manifest_path.parent),
        (manifest["remote_project"] + "/logs", manifest_path.parent / "logs"),
        (manifest["remote_project"] + "/provenance", manifest_path.parent / "environments"),
    ):
        destination.mkdir(exist_ok=True)
        log_filters = []
        if destination.name == "logs":
            # Copy only tracked job logs; previous monitors duplicated every project log.
            log_filters = [f"--include=*_{job_id}.log" for job_id in ids] + ["--exclude=*"]
        subprocess.run(
            [
                "rsync",
                "-az",
                "--timeout=60",
                "--exclude=/manifest.json",
                *log_filters,
                "-e",
                shlex.join(SSH),
                f"hoffman2:{source}/",
                str(destination) + "/",
            ],
            check=True,
            capture_output=True,
            text=True,
            timeout=180,
        )
    print(
        snapshot["checked_at"],
        "; ".join(f"{name}={state['state']}" for name, state in snapshot["jobs"].items()),
        flush=True,
    )
    try:
        refresh_prediction_review(manifest_path, manifest)
    except (subprocess.SubprocessError, OSError) as error:
        print(
            f"Review refresh failed; scheduler monitoring continues: {error}",
            file=sys.stderr,
            flush=True,
        )
    return done


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--watch", action="store_true")
    parser.add_argument("--interval", type=int, default=120)
    args = parser.parse_args()
    if args.interval < 60:
        parser.error("Poll no more than once a minute")
    manifest_path = args.manifest.resolve()
    while True:
        try:
            done = check(manifest_path)
            if done:
                manifest = json.loads(manifest_path.read_text())
                if manifest.get("workflow") == "cleavage_protector_v1":
                    subprocess.run(
                        [
                            sys.executable,
                            str(manifest_path.parent / "workflow/protector_round.py"),
                            "report",
                            str(manifest_path.parent),
                        ],
                        check=True,
                    )
                    return
                if manifest.get("experiments"):
                    subprocess.run(
                        [
                            sys.executable,
                            str(manifest_path.parent / "workflow/bindcraft_round.py"),
                            "report",
                            str(manifest_path.parent),
                        ],
                        check=True,
                    )
                    return
                subprocess.run(
                    [
                        sys.executable,
                        str(ROOT / "scripts/plot_pilot.py"),
                        str(manifest_path.parent),
                    ],
                    check=True,
                )
            if done or not args.watch:
                return
        except (subprocess.SubprocessError, OSError, ValueError) as error:
            print(f"Monitor error: {error}", file=sys.stderr, flush=True)
            if not args.watch:
                raise
        time.sleep(args.interval)


if __name__ == "__main__":
    main()
