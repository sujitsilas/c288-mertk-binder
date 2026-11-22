#!/usr/bin/env python3
"""Execute and report a dedicated stalk-protector pilot with explicit uncertainty."""

import argparse
import json
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from ig_followup import report_af3, run_af3, write

from mertk.campaign_metrics import summarize
from mertk.cleavage_protector import evaluate_pose
from mertk.protector_screen import rank_key, saved_poses
from mertk.provenance import file_hash


def preflight(suite, gpu):
    record = json.loads((suite / "manifest.json").read_text())
    py = Path(record["remote_project"]) / "envs/bindcraft-wheels/bin/python"
    for key in record["experiments"]:
        inp = suite / "campaigns" / key / "inputs"
        command = [
            str(py),
            str(suite / "workflow/check_glycan_runtime.py"),
            "--target",
            str(inp / "target.pdb"),
            "--constraints",
            str(inp / "glycan_constraints.json"),
            "--hotspot",
            "A490",
            "--output",
            str(suite / f"{'gpu' if gpu else 'cpu'}_preflight_{key}.json"),
        ]
        if gpu:
            command += [
                "--gpu-forward",
                "--params",
                record["remote_project"] + "/external/BindCraft",
            ]
        subprocess.run(
            command, env=dict(os.environ, JAX_PLATFORMS="cuda" if gpu else "cpu"), check=True
        )


def check_preflight(suite, key):
    inp = suite / "campaigns" / key / "inputs"
    q = json.loads((suite / f"gpu_preflight_{key}.json").read_text())
    if (
        not q.get("passed")
        or not q.get("full_forward_backward")
        or q["constraints_sha256"] != file_hash(inp / "glycan_constraints.json")
        or q["target_sha256"] != file_hash(inp / "target.pdb")
        or "membrane_volume" not in q.get("model_losses", {})
    ):
        raise ValueError("Matched full GPU glycan/membrane preflight required")
    if not json.loads((suite / "workflow/positive_control_evidence.json").read_text())["passed"]:
        raise ValueError("Positive-control evidence missing")
    print("PROTECTOR_GPU_PREFLIGHT_VERIFIED", key, flush=True)


def quality(suite, key):
    campaign = suite / "campaigns" / key
    metrics = summarize(campaign)
    rows = []
    exploratory = (
        json.loads((suite / "manifest.json").read_text()).get("protector_screen") == "exploratory"
    )
    for path, meta in saved_poses(campaign / "bindcraft/designs"):
        if not exploratory and path.parent.name != "Accepted":
            continue
        try:
            if "error" in meta:
                raise ValueError(meta["error"])
            row = evaluate_pose(path, campaign, suite)
            if meta["expected_sequence"] and row["binder_sequence"] != meta["expected_sequence"]:
                raise ValueError("Statistics and coordinate sequences differ")
            row.update(meta, coordinate_file=str(path.relative_to(campaign)))
            if path.parent.name == "Accepted":
                row["mean_af2_iptm"] = meta["iptm"]
        except (ValueError, KeyError, IndexError) as error:
            row = dict(path=str(path), preliminary_model_pass=False, error=str(error))
        rows.append(row)
    result = dict(
        status="Exploratory model/context screen; no native-binding or cleavage-protection claim",
        records=rows,
        accepted_variants=sum(Path(r["path"]).parent.name == "Accepted" for r in rows),
        reviewable_poses=sum("error" not in r for r in rows),
        preliminary_model_passes=sum(r["preliminary_model_pass"] for r in rows),
        qualified=0,
        cleavage_protection_predicted=False,
    )
    write(campaign / "metrics.json", metrics)
    write(campaign / "protector_quality.json", result)
    print(
        json.dumps(
            {k: result[k] for k in ["accepted_variants", "preliminary_model_passes", "qualified"]}
        ),
        flush=True,
    )


def select_parents(suite, *, manage_scheduler=True):
    record = json.loads((suite / "manifest.json").read_text())
    exploratory = record.get("protector_screen") == "exploratory"
    selected, reviewed = [], []
    for key in record["experiments"]:
        campaign = suite / "campaigns" / key
        path = campaign / "protector_quality.json"
        if not path.exists():
            reviewed.append({"experiment": key, "error": "Missing quality report"})
            continue
        q = json.loads(path.read_text())
        reviewed.append(
            {
                "experiment": key,
                "accepted": q["accepted_variants"],
                "model_passes": q["preliminary_model_passes"],
            }
        )
        for row in q["records"]:
            if "error" in row:
                continue
            if exploratory:
                parent = campaign / row["coordinate_file"]
                order = rank_key(row)
            elif row.get("preliminary_model_pass"):
                parent = campaign / "bindcraft/designs/Accepted" / Path(row["path"]).name
                order = (-row["mean_af2_iptm"],)
            else:
                continue
            if file_hash(parent) != row["sha256"]:
                raise ValueError("Selected parent changed after screening")
            selected.append((order, key, parent, row))
    selected.sort(key=lambda x: (x[0], x[1], x[2].name))
    # One parent total; AF3 budget is reserved even while its job is held.
    cases = []
    if selected and not any("error" in r for r in reviewed):
        _, key, parent, row = selected[0]
        name = "mertk_stalk_protector_parent"
        data = dict(
            name=name,
            dialect="alphafold3",
            version=1,
            modelSeeds=[101, 202, 303],
            sequences=[
                {
                    "protein": {
                        "id": c,
                        "sequence": s,
                        "unpairedMsa": "",
                        "pairedMsa": "",
                        "templates": [],
                    }
                }
                for c, s in [
                    ("A", record["canonical_target_sequence"]),
                    ("B", row["binder_sequence"]),
                ]
            ],
        )
        path = suite / "assessment/inputs" / f"{name}.json"
        write(path, data)
        cases.append(
            dict(
                id=name,
                name=name,
                input=f"inputs/{name}.json",
                input_sha256=file_hash(path),
                sequence=row["binder_sequence"],
                epitope="Stalk_P490S491",
                experiment=key,
                parent=str(parent.relative_to(suite)),
                parent_sha256=file_hash(parent),
                parent_stage=row.get("stage", "Sequence refold"),
                status="Exploratory model-screened parent; unvalidated protector",
            )
        )
    manifest = dict(
        created_at=datetime.now(timezone.utc).isoformat(),
        cases=cases,
        reviewed=reviewed,
        protocol=record["af3_protocol"],
        target={"sequence": record["canonical_target_sequence"], "native_range": [93, 526]},
        comparator_note="Target includes stalk/TM. Earlier Ig controls used A93–A484; scores are not a matched comparison.",
    )
    write(suite / "assessment/manifest.json", manifest)
    if not cases:
        write(
            suite / "assessment/execution_skipped.json",
            {
                "reason": "No eligible model-screened stalk parent; no AF3 GPU allocation required",
                "reviewed": reviewed,
            },
        )
    if manage_scheduler:
        job = record["jobs"]["af3_assessment"]["job_id"]
        before = subprocess.check_output(["scontrol", "show", "job", job, "-o"], text=True)
        if (
            "JobState=PENDING " not in before
            or record["campaign_id"] + "_af3_assessment" not in before
        ):
            raise ValueError("Expected the held assessment job for this suite")
        command = ["scontrol", "release", job] if cases else ["scancel", "--state=PENDING", job]
        subprocess.run(command, check=True)
        write(
            suite / "assessment/scheduler_decision.json",
            {
                "at": datetime.now(timezone.utc).isoformat(),
                "before": before,
                "command": command,
                "eligible_parents": len(cases),
            },
        )
    print(f"PROTECTOR_AF3_PARENTS {len(cases)}", flush=True)


def report(suite):
    record = json.loads((suite / "manifest.json").read_text())
    report_af3(suite)
    summary = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "qualified_protectors": 0,
        "cleavage_protection_tested": False,
        "experiments": {},
        "required_experiments": [
            "MERTK binding",
            "Reduced soluble MERTK with retained surface receptor",
            "Preserved GAS6/PROS1 function and efferocytosis",
            "ADAM17 activity/selectivity control",
        ],
    }
    for key in record["experiments"]:
        q = suite / "campaigns" / key / "protector_quality.json"
        summary["experiments"][key] = (
            json.loads(q.read_text()) if q.exists() else {"status": "Quality pending or absent"}
        )
    write(suite / "protector_summary.json", summary)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "stage",
        choices=["validate", "preflight", "check-preflight", "quality", "select", "af3", "report"],
    )
    parser.add_argument("suite", type=Path)
    parser.add_argument("--experiment", default="")
    args = parser.parse_args()
    if args.stage in {"validate", "preflight"}:
        preflight(args.suite, args.stage == "preflight")
    elif args.stage == "check-preflight":
        check_preflight(args.suite, args.experiment)
    elif args.stage == "quality":
        quality(args.suite, args.experiment)
    elif args.stage == "select":
        select_parents(args.suite)
    elif args.stage == "af3":
        os.environ["XLA_FLAGS"] = "--xla_gpu_enable_triton_gemm=false"
        record = json.loads((args.suite / "manifest.json").read_text())
        subprocess.run(
            [
                record["remote_project"] + "/envs/alphafold3/bin/python",
                "-c",
                "import jax; print(jax.devices()); assert any(d.platform == 'gpu' for d in jax.devices())",
            ],
            check=True,
        )
        run_af3(args.suite, "parents")
    else:
        report(args.suite)


if __name__ == "__main__":
    main()
