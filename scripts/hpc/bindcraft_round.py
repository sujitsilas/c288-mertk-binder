#!/usr/bin/env python3
"""Validate frozen inputs, gate the control, and report real BindCraft outputs."""

import argparse
import json
import subprocess
from pathlib import Path

from check_course_outputs import analyze, collect

from mertk.campaign_metrics import summarize
from mertk.provenance import file_hash


def validate(suite, *, remote=False):
    record = json.loads((suite / "manifest.json").read_text())
    for relative, expected in record["input_sha256"].items():
        if file_hash(suite / relative) != expected:
            raise ValueError(f"Frozen input/workflow changed: {relative}")
    expected_filters = record["bindcraft_default_filter_sha256"]
    if (
        record.get("workflow") == "cleavage_protector_v1"
        and record.get("protector_screen") == "exploratory"
    ):
        expected_filters = record["protector_filter_sha256"]
    for key in record["experiments"]:
        filters = suite / "campaigns" / key / "inputs/bindcraft/filters.json"
        if file_hash(filters) != expected_filters:
            raise ValueError("Frozen BindCraft filters changed")
    if remote:
        project = Path(record["remote_project"])
        engine = project / "external/BindCraft"
        revision = subprocess.check_output(
            ["git", "-C", str(engine), "rev-parse", "HEAD"], text=True
        ).strip()
        if revision != record["bindcraft_revision"]:
            raise ValueError("Remote engine revision changed")
        # Compiled platform-specific binaries may differ; Python sources must be pristine.
        subprocess.run(
            ["git", "-C", str(engine), "diff", "--exit-code", "HEAD", "--", "*.py"], check=True
        )
        ready = json.loads((project / "provenance/bindcraft-wheels-ready.json").read_text())
        if ready.get("gpu_validated") is not True:
            raise ValueError("BindCraft GPU runtime not validated")
    print("FROZEN_ROUND_INPUTS_VERIFIED", flush=True)
    return record


def report_experiment(campaign, control=False):
    manifest = json.loads((campaign / "manifest.json").read_text())
    metrics = summarize(campaign)
    paths = sorted((campaign / "bindcraft/designs/Accepted").glob("*.pdb"))
    if control:
        rows = []
        for path in paths:
            try:
                row = analyze(path, campaign / "inputs/target.pdb")
            except (ValueError, KeyError, IndexError) as error:
                row = {"path": str(path), "preliminary_geometry_pass": False, "error": str(error)}
            rows.append(row)
        passed = (
            metrics["accepted_variants"] > 0
            and not metrics["accepted_evidence_mismatch"]
            and any(r["preliminary_geometry_pass"] for r in rows)
        )
        result = {
            "passed": passed,
            "metrics": metrics,
            "geometry": rows,
            "interpretation": "Diagnostic positive control only; zero hits is inconclusive, not proof of a broken stack.",
        }
        (campaign / "control_gate.json").write_text(json.dumps(result, indent=2) + "\n")
    else:
        rows = collect(
            paths,
            campaign / "inputs/target.pdb",
            campaign / "inputs/37NT.cif",
            hotspots=manifest["hotspots"],
        )
        result = {
            "status": "Accepted-only preliminary native-context screening; manual sequon and AF3 review still required",
            "arms": {
                "bindcraft": {
                    "records": rows,
                    "files": len(rows),
                    "preliminary_passes": sum(r["preliminary_geometry_pass"] for r in rows),
                }
            },
        }
        (campaign / "native_quality.json").write_text(json.dumps(result, indent=2) + "\n")
        passed = True
    (campaign / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    print(
        json.dumps(
            {
                k: metrics[k]
                for k in [
                    "campaign_id",
                    "attempts_started",
                    "attempts_completed",
                    "accepted_variants",
                    "accepted_trajectories",
                ]
            }
        ),
        flush=True,
    )
    return passed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["validate", "control-gate", "report"])
    parser.add_argument("suite", type=Path)
    parser.add_argument("--experiment")
    parser.add_argument("--remote", action="store_true")
    args = parser.parse_args()
    record = validate(args.suite, remote=args.remote)
    if args.stage == "control-gate":
        if not report_experiment(args.suite / "campaigns/pdl1_control", control=True):
            raise SystemExit(42)
    elif args.stage == "report":
        keys = [args.experiment] if args.experiment else list(record["experiments"])
        for key in keys:
            if key not in record["experiments"]:
                raise ValueError("Unknown experiment")
            report_experiment(args.suite / "campaigns" / key, control=key == "pdl1_control")


if __name__ == "__main__":
    main()
