#!/usr/bin/env python3
"""Run a prepared AF3 bundle in a separate environment, with resumable provenance."""

import argparse
import json
import shutil
import subprocess
import time
from pathlib import Path

from mertk.af3 import load_jobs
from mertk.structural import prediction_files, validate_processed_input
from mertk.training import file_hash


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--af3-python", type=Path, required=True)
    parser.add_argument("--model-dir", type=Path, required=True)
    parser.add_argument("--af3-root", type=Path, default=Path("external/alphafold3"))
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    manifest = load_jobs(args.manifest)
    bundle, af3_root = args.manifest.resolve().parent, args.af3_root.resolve()
    predictions = bundle / "predictions"
    protocol = manifest["protocol"]
    command = [
        str(args.af3_python.resolve()),
        str(af3_root / "run_alphafold.py"),
        f"--input_dir={bundle / 'pending_inputs'}",
        f"--output_dir={predictions}",
        f"--model_dir={args.model_dir.resolve()}",
        "--run_data_pipeline=false",
        "--run_inference=true",
        f"--num_diffusion_samples={protocol['samples_per_seed']}",
        f"--num_recycles={protocol['num_recycles']}",
        f"--flash_attention_implementation={protocol['flash_attention_implementation']}",
        "--compress_large_output_files=false",
    ]
    if args.resume:
        command.append("--force_output_dir=true")
    if not args.execute:
        print(
            json.dumps(
                {"command": command, "planned_candidates": len(manifest["jobs"]), "execute": False},
                indent=2,
            )
        )
        return
    revision = subprocess.check_output(
        ["git", "-C", str(af3_root), "rev-parse", "HEAD"], text=True
    ).strip()
    dirty = subprocess.check_output(
        ["git", "-C", str(af3_root), "status", "--porcelain", "--untracked-files=no"], text=True
    )
    if revision != manifest["engine"]["revision"] or dirty:
        raise ValueError("AF3 checkout differs from the source pinned in this plan")
    weights = sorted(path for path in args.model_dir.resolve().glob("*.bin*") if path.is_file())
    if not weights:
        raise FileNotFoundError("AF3 model directory has no weights")
    environment = json.loads(
        subprocess.check_output(
            [
                str(args.af3_python.resolve()),
                "-c",
                "import importlib.metadata as m,json,platform; print(json.dumps({'python':platform.python_version(),'packages':sorted((d.metadata['Name'],d.version) for d in m.distributions())}))",
            ],
            text=True,
        )
    )
    engine = {
        "revision": revision,
        "environment": environment,
        "weights": {
            str(path.relative_to(args.model_dir.resolve())): file_hash(path) for path in weights
        },
    }
    record_path = predictions / "execution.json"
    record = {"plan_id": manifest["plan_id"], "engine": engine, "attempts": []}
    if record_path.exists():
        previous = json.loads(record_path.read_text())
        if (
            not args.resume
            or previous["plan_id"] != record["plan_id"]
            or previous["engine"] != engine
        ):
            raise ValueError("Existing AF3 run requires --resume with identical plan and engine")
        record = previous
    elif predictions.exists() and any(predictions.iterdir()):
        raise ValueError("Untracked predictions already exist; choose a fresh bundle")
    pending = []
    for job in manifest["jobs"]:
        if args.resume:
            try:
                prediction_files(predictions, job, protocol)
                validate_processed_input(
                    predictions / job["name"] / f"{job['name']}_data.json",
                    json.loads((bundle / job["input"]).read_text()),
                )
                continue
            except FileNotFoundError:
                pass
        pending.append(job)
    if not pending:
        print("All planned AF3 outputs are present; run scripts/score_af3.py.")
        return
    pending_dir = bundle / "pending_inputs"
    pending_dir.mkdir(exist_ok=True)
    # Only remove our previous staging copies, never prediction results.
    for path in pending_dir.glob("*.json"):
        path.unlink()
    for job in pending:
        shutil.copy2(bundle / job["input"], pending_dir / Path(job["input"]).name)
    predictions.mkdir(exist_ok=True)
    attempt = {
        "candidates": [job["name"] for job in pending],
        "command": command,
        "status": "running",
    }
    record["attempts"].append(attempt)
    record_path.write_text(json.dumps(record, indent=2) + "\n")
    started = time.perf_counter()
    result = subprocess.run(command)
    attempt.update(
        status="completed" if result.returncode == 0 else "failed",
        returncode=result.returncode,
        seconds=time.perf_counter() - started,
    )
    record_path.write_text(json.dumps(record, indent=2) + "\n")
    result.check_returncode()


if __name__ == "__main__":
    main()
