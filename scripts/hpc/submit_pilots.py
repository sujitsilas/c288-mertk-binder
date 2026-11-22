#!/usr/bin/env python3
"""Submit one setup and dependent GPU pilot per arm, recording IDs immediately."""

import argparse
import json
import shlex
import subprocess
from pathlib import Path

from mertk.target_policy import check_manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--account", required=True, help="User-selected Slurm account")
    args = parser.parse_args()
    record = check_manifest(args.manifest, generation=True)
    if record["jobs"]:
        raise SystemExit("Jobs already recorded; inspect them before starting another campaign.")
    ssh = [
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
    project, campaign = record["remote_project"], record["remote_campaign"]
    for tool in ("bindcraft", "rfd4_proteina", "alphafold3"):
        dependency = None
        # AF3 evaluates real candidates later; no empty inference job is queued.
        stages = ("setup",) if tool == "alphafold3" else ("setup", "generate")
        for stage in stages:
            name = f"{tool}_{stage}"
            cmd = [
                "sbatch",
                "--parsable",
                f"--account={args.account}",
                "--nodes=1",
                "--cpus-per-task=8",
                f"--job-name=mertk_{name}",
                f"--output={project}/logs/%x_%j.log",
                f"--chdir={project}",
            ]
            if stage == "setup":
                cmd += [
                    "--time=04:00:00",
                    "--mem=32G",
                    f"{project}/scripts/hpc/setup_runtime.sh",
                    tool,
                    project,
                ]
            else:
                cmd += [
                    "--time=02:00:00",
                    "--mem=64G",
                    "--gres=gpu:1",
                    "--constraint=A100",
                    f"--dependency=afterok:{dependency}",
                    "--kill-on-invalid-dep=yes",
                    f"{project}/scripts/hpc/run_pilot.sh",
                    tool,
                    project,
                    campaign,
                ]
            result = subprocess.run(
                ssh + [shlex.join(cmd)], check=True, text=True, stdout=subprocess.PIPE
            )
            job_id = result.stdout.strip().split(";")[0]
            if not job_id.isdecimal():
                raise SystemExit(f"Unexpected submission response: {result.stdout!r}")
            record["jobs"][name] = {"job_id": job_id, "command": cmd}
            args.manifest.write_text(json.dumps(record, indent=2) + "\n")
            print(f"{name}: {job_id}", flush=True)
            dependency = job_id


if __name__ == "__main__":
    main()
