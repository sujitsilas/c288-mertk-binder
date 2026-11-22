#!/usr/bin/env python3
"""Matched AF3 pilot and monovalent protector triage; never infer biological activity."""

import argparse
import csv
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from mertk.provenance import file_hash
from mertk.structural import prediction_files, read_prediction, validate_processed_input


def write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n")


def prepare_candidates(suite):
    record = json.loads((suite / "manifest.json").read_text())
    controls = json.loads((suite / "comparisons/manifest.json").read_text())
    selected, reviewed = [], []
    for key in record["experiments"]:
        campaign = suite / "campaigns" / key
        quality_path = campaign / "native_quality.json"
        if not quality_path.exists():
            reviewed.append({"epitope": key, "status": "Missing native quality report"})
            continue
        quality = json.loads(quality_path.read_text())["arms"]["bindcraft"]["records"]
        stats_path = campaign / "bindcraft/designs/final_design_stats.csv"
        scores = {}
        if stats_path.exists():
            with stats_path.open() as handle:
                scores = {r["Design"]: float(r["Average_i_pTM"]) for r in csv.DictReader(handle)}
        good = [q for q in quality if q.get("preliminary_geometry_pass")]
        good.sort(
            key=lambda q: (-scores.get(Path(q["path"]).stem.rsplit("_model", 1)[0], -1), q["path"])
        )
        reviewed.append({"epitope": key, "accepted": len(quality), "native_passes": len(good)})
        if not good:
            continue
        q = good[0]
        parent = Path(q["path"])
        # Local reports may contain local absolute paths; use only the known campaign basename.
        parent = campaign / "bindcraft/designs/Accepted" / parent.name
        if file_hash(parent) != q["sha256"]:
            raise ValueError("Selected parent changed after native screening")
        name = f"mertk_{key.lower()}_parent"
        input_data = {
            "name": name,
            "dialect": "alphafold3",
            "version": 1,
            "modelSeeds": controls["protocol"]["model_seeds"],
            "sequences": [
                {
                    "protein": {
                        "id": chain,
                        "sequence": seq,
                        "unpairedMsa": "",
                        "pairedMsa": "",
                        "templates": [],
                    }
                }
                for chain, seq in [
                    ("A", controls["target"]["sequence"]),
                    ("B", q["binder_sequence"]),
                ]
            ],
        }
        path = suite / "assessment/inputs" / f"{name}.json"
        write(path, input_data)
        selected.append(
            {
                "id": name,
                "name": name,
                "input": str(path.relative_to(suite / "assessment")),
                "input_sha256": file_hash(path),
                "sequence": q["binder_sequence"],
                "epitope": key,
                "parent": str(parent.relative_to(suite)),
                "parent_sha256": q["sha256"],
                "status": "Provisional computational triage only; not a native-binding nomination",
            }
        )
    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "cases": selected,
        "protocol": controls["protocol"],
        "target": controls["target"],
        "reviewed": reviewed,
    }
    write(suite / "assessment/manifest.json", manifest)
    write(
        suite / "protector/triage.json",
        {
            "status": "Monovalent parents await experimental shedding screen"
            if selected
            else "No native-context-passing parents available",
            "shared_library": True,
            "parents": selected,
            "reviewed": reviewed,
            "cleavage_protection_predicted": False,
            "reason": "Ig1/Ig2 binding is distal to the reported P490-S491 bond; AF3 cannot establish inhibition of proteolysis. 37NT ends at G484.",
            "format": "Monovalent; separate from receptor-crosslinking formats",
            "required_measurements": [
                "Soluble MERTK release",
                "Intact surface and total MERTK",
                "Ligand-supported internalized apoptotic cargo",
                "Assay interference and receptor depletion controls",
            ],
            "dimerization_formats": [
                "Two identical Ig1 arms",
                "Two identical Ig2 arms",
                "Ig1-Ig2 tandem tested for one-receptor versus two-receptor occupancy",
            ],
            "architecture_status": "No fused sequence or active geometry nominated",
        },
    )
    print(
        f"Prepared {len(selected)} provisional parents for matched AF3; protector triage recorded",
        flush=True,
    )


def run_af3(suite, case_id):
    record = json.loads((suite / "manifest.json").read_text())
    project = Path(record["remote_project"])
    bundle = suite / ("assessment" if case_id == "parents" else "comparisons")
    manifest = json.loads((bundle / "manifest.json").read_text())
    cases = (
        manifest["cases"]
        if case_id == "parents"
        else [c for c in manifest["cases"] if c["id"] == case_id]
    )
    if not cases:
        write(
            bundle / "execution_skipped.json",
            {"reason": "No native-context-passing parent; GPU inference skipped"},
        )
        print("NO_QUALIFIED_PARENT_AF3_SKIPPED", flush=True)
        return
    tools = json.loads((suite / "workflow/external_tools.json").read_text())
    engine = project / tools["alphafold3"]["directory"]
    revision = subprocess.check_output(
        ["git", "-C", str(engine), "rev-parse", "HEAD"], text=True
    ).strip()
    if revision != tools["alphafold3"]["revision"]:
        raise ValueError("AF3 source revision differs from frozen plan")
    subprocess.run(
        ["git", "-C", str(engine), "diff", "--exit-code", "HEAD", "--", "*.py"], check=True
    )
    weights = project / "weights/alphafold3/af3.bin.zst"
    if file_hash(weights) != record["af3_weights_sha256"]:
        raise ValueError("AF3 weights differ from verified download")
    protocol = manifest["protocol"]
    for case in cases:
        inp = bundle / case["input"]
        if file_hash(inp) != case["input_sha256"]:
            raise ValueError("AF3 input checksum mismatch")
        expected = json.loads(inp.read_text())
        if expected["modelSeeds"] != [101, 202, 303]:
            raise ValueError("Matched three-seed protocol required")
        output = bundle / "predictions" / case["id"]
        if output.exists():
            raise ValueError("Preserve previous AF3 attempt; do not silently overwrite")
        command = [
            str(project / "envs/alphafold3/bin/python"),
            str(engine / "run_alphafold.py"),
            f"--json_path={inp}",
            f"--output_dir={bundle / 'predictions'}",
            f"--model_dir={weights.parent}",
            "--run_data_pipeline=false",
            "--run_inference=true",
            "--num_diffusion_samples=1",
            "--num_recycles=10",
            "--flash_attention_implementation=triton",
            "--compress_large_output_files=false",
        ]
        execution = {
            "case": case["id"],
            "command": command,
            "input_sha256": case["input_sha256"],
            "engine_revision": revision,
            "weights_sha256": record["af3_weights_sha256"],
            "protocol": protocol,
            "started_at": datetime.now(timezone.utc).isoformat(),
            "status": "running",
        }
        path = bundle / "execution" / f"{case['id']}.json"
        write(path, execution)
        completed = subprocess.run(command)
        execution.update(
            returncode=completed.returncode,
            status="completed" if completed.returncode == 0 else "failed",
        )
        write(path, execution)
        completed.check_returncode()


def report_af3(suite):
    for name in ("comparisons", "assessment"):
        bundle = suite / name
        if not (bundle / "manifest.json").exists():
            continue
        manifest = json.loads((bundle / "manifest.json").read_text())
        errors = []
        for case in manifest["cases"]:
            case["predictions"] = []
            try:
                inp = bundle / case["input"]
                if file_hash(inp) != case["input_sha256"]:
                    raise ValueError("Input checksum mismatch")
                expected = json.loads(inp.read_text())
                directory = bundle / "predictions" / case["id"]
                validate_processed_input(directory / f"{case['id']}_data.json", expected)
                samples = prediction_files(
                    bundle / "predictions", {"name": case["id"]}, manifest["protocol"]
                )
                predictions = []
                for sample_name, paths in samples:
                    read_prediction(
                        paths["model"], manifest["target"]["sequence"], case["sequence"]
                    )
                    summary = json.loads(paths["summary"].read_text())
                    confidences = json.loads(paths["confidence"].read_text())
                    chain = np.asarray(confidences["atom_chain_ids"])
                    plddt = np.asarray(confidences["atom_plddts"])
                    pair = np.asarray(summary["chain_pair_iptm"], dtype=float)
                    if pair.shape != (2, 2) or not np.isfinite(pair).all():
                        raise ValueError("Invalid two-chain ipTM matrix")
                    predictions.append(
                        {
                            "pose": str(paths["model"].relative_to(bundle)),
                            "pose_sha256": file_hash(paths["model"]),
                            "input_sha256": case["input_sha256"],
                            "seed": int(sample_name.split("_")[0].split("-")[1]),
                            "iptm": float(pair[0, 1]),
                            "partner_plddt": float(plddt[chain == "B"].mean()),
                            "ranking_score": summary["ranking_score"],
                            "has_clash": summary.get("has_clash"),
                        }
                    )
                case.update(
                    predictions=predictions,
                    status="All planned seeds validated; prediction is not binding evidence",
                )
            except (FileNotFoundError, ValueError, KeyError) as error:
                case["status"] = f"Incomplete or invalid prediction: {error}"
                errors.append({"case": case["id"], "error": str(error)})
        manifest.update(reported_at=datetime.now(timezone.utc).isoformat(), errors=errors)
        write(bundle / "results_manifest.json", manifest)
        print(
            f"{name}: {sum(bool(c['predictions']) for c in manifest['cases'])}/{len(manifest['cases'])} complete cases",
            flush=True,
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["prepare", "af3", "report"])
    parser.add_argument("suite", type=Path)
    parser.add_argument("--case")
    args = parser.parse_args()
    if args.stage == "prepare":
        prepare_candidates(args.suite)
    elif args.stage == "af3":
        run_af3(args.suite, args.case)
    else:
        report_af3(args.suite)


if __name__ == "__main__":
    main()
