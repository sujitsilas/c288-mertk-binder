import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from Bio.PDB import MMCIFIO, Atom, Chain, Model, Residue, Structure
from Bio.SeqUtils import seq3

from mertk.af3 import load_jobs, prepare_jobs
from mertk.structural import geometry_metrics, rigid_fit, score_bundle, summarize_structures
from mertk.training import file_hash


@pytest.fixture
def candidate_fixture(tmp_path):
    n, m = 8, 12
    t = np.arange(m, dtype=float)
    target_ca = np.column_stack((1.3 * t, 2 * np.sin(t), 2 * np.cos(t)))
    binder_ca = target_ca[:n] + [0, 4, 0]
    mask = np.array([False] + [True] * (n - 1))
    reference = {
        "binder_sequence": np.array("A" * n),
        "target_sequence": np.array("G" * m),
        "binder_backbone": np.stack((binder_ca + [-1, 1, 0], binder_ca, binder_ca + [1, 0, 0]), 1),
        "target_backbone": np.stack((target_ca + [-1, 1, 0], target_ca, target_ca + [1, 0, 0]), 1),
        "design_mask": mask,
    }
    source = tmp_path / "reference.npz"
    np.savez_compressed(source, **reference)
    base = {
        "complex_id": "fixture",
        "component": "family",
        "split": "test",
        "reference_npz": str(source),
        "reference_sha256": file_hash(source),
        "mask_sha256": hashlib.sha256(mask.astype(np.uint8).tobytes()).hexdigest(),
        "dataset_identity": "fixture-only",
        "target_sequence": "G" * m,
        "binder_sequence": "A" * n,
        "design_positions_1based": list(range(2, n + 1)),
        "reference_kind": "native",
        "smoke_test": False,
        "control": False,
    }
    arms = ["frozen_esm3", "lora", "combined", "proteinmpnn"]
    candidates = [
        {**base, "candidate_id": f"{arm}:42:{i}", "arm": arm, "seed": 42, "sample_index": i}
        for arm in arms
        for i in range(2)
    ]
    candidates.append(
        {
            **base,
            "candidate_id": "original",
            "arm": "native_original",
            "seed": None,
            "sample_index": 0,
            "control": True,
        }
    )
    path = tmp_path / "candidates.jsonl"
    path.write_text("".join(json.dumps(row) + "\n" for row in candidates))
    protocol = json.loads(Path("configs/af3_evaluation.json").read_text())
    protocol.update(target_msa="single_sequence", samples_per_seed=2)
    return path, protocol, reference, candidates


def write_prediction(path, target_sequence, binder_sequence, target_ca, binder_ca):
    structure = Structure.Structure("fixture")
    model = Model.Model(0)
    structure.add(model)
    for chain_id, sequence, coords in [
        ("A", target_sequence, target_ca),
        ("B", binder_sequence, binder_ca),
    ]:
        chain = Chain.Chain(chain_id)
        model.add(chain)
        for i, (aa, coord) in enumerate(zip(sequence, coords, strict=True), 1):
            residue = Residue.Residue((" ", i, " "), seq3(aa).upper(), " ")
            residue.add(
                Atom.Atom("CA", coord.astype(np.float32), 90, 1, " ", " CA ", i, element="C")
            )
            chain.add(residue)
    writer = MMCIFIO()
    writer.set_structure(structure)
    writer.save(str(path))


def fake_outputs(bundle, plan, reference, bad_top_pose=False):
    predictions = bundle / "predictions"
    for job in plan["jobs"]:
        directory = predictions / job["name"]
        directory.mkdir(parents=True)
        (directory / f"{job['name']}_data.json").write_text((bundle / job["input"]).read_text())
        candidate = job["candidate"]
        target_ca, binder_ca = (
            reference["target_backbone"][:, 1],
            reference["binder_backbone"][:, 1],
        )
        chains = ["A"] * len(target_ca) + ["B"] * len(binder_ca)
        with (directory / f"{job['name']}_ranking_scores.csv").open("w") as handle:
            writer = csv.writer(handle)
            writer.writerow(["seed", "sample", "ranking_score"])
            writer.writerows(
                (seed, sample, 0.9 - 0.1 * sample)
                for seed in plan["protocol"]["model_seeds"]
                for sample in range(plan["protocol"]["samples_per_seed"])
            )
        for seed in plan["protocol"]["model_seeds"]:
            for sample in range(plan["protocol"]["samples_per_seed"]):
                sub = directory / f"seed-{seed}_sample-{sample}"
                sub.mkdir()
                summary = {
                    "ranking_score": 0.9 - 0.1 * sample,
                    "chain_ids": ["B", "A"],
                    "chain_pair_iptm": [[0.8, 0.81], [0.9, 0.85]],
                    "has_clash": 0.0,
                }
                (sub / "fixture_summary_confidences.json").write_text(json.dumps(summary))
                confidence = {
                    "token_chain_ids": chains,
                    "token_res_ids": list(range(1, len(target_ca) + 1))
                    + list(range(1, len(binder_ca) + 1)),
                    "atom_chain_ids": chains,
                    "atom_plddts": [90] * len(chains),
                    "pae": np.full((len(chains), len(chains)), 2.0).tolist(),
                }
                (sub / "fixture_confidences.json").write_text(json.dumps(confidence))
                translated = binder_ca + ([0, 0, 6] if sample == 0 and bad_top_pose else [0, 0, 0])
                write_prediction(
                    sub / "fixture_model.cif",
                    candidate["target_sequence"],
                    candidate["binder_sequence"],
                    target_ca,
                    translated,
                )
    return predictions


def test_preparation_uses_shared_msas_and_no_templates(candidate_fixture, tmp_path):
    path, protocol, _, candidates = candidate_fixture
    protocol["target_msa"] = "provided"
    msa = tmp_path / "target.a3m"
    msa.write_text(">target\n" + candidates[0]["target_sequence"] + "\n")
    key = hashlib.sha256(candidates[0]["target_sequence"].encode()).hexdigest()
    plan = prepare_jobs(
        path, protocol, tmp_path / "bundle", {key: str(msa)}, samples_per_design_seed=1
    )
    assert len(plan["jobs"]) == 5  # Four matched redesigns plus the original control.
    target_inputs = []
    for job in plan["jobs"]:
        data = json.loads((tmp_path / "bundle" / job["input"]).read_text())
        target, binder = [row["protein"] for row in data["sequences"]]
        assert target["templates"] == binder["templates"] == []
        assert binder["unpairedMsa"] == binder["pairedMsa"] == ""
        target_inputs.append(target)
    assert all(target == target_inputs[0] for target in target_inputs)
    assert load_jobs(tmp_path / "bundle/manifest.json")["plan_id"] == plan["plan_id"]
    (tmp_path / "bundle" / plan["jobs"][0]["input"]).write_text("{}")
    with pytest.raises(ValueError, match="input changed"):
        load_jobs(tmp_path / "bundle/manifest.json")


def test_af3_launcher_can_plan_without_weights_or_running_inference(candidate_fixture, tmp_path):
    path, protocol, _, _ = candidate_fixture
    bundle = tmp_path / "bundle"
    prepare_jobs(path, protocol, bundle)
    result = json.loads(
        subprocess.check_output(
            [
                sys.executable,
                "scripts/run_af3.py",
                "--manifest",
                str(bundle / "manifest.json"),
                "--af3-python",
                sys.executable,
                "--model-dir",
                str(tmp_path / "absent_weights"),
            ],
            text=True,
        )
    )
    assert result["execute"] is False
    assert result["planned_candidates"] == 9
    assert "--run_data_pipeline=false" in result["command"]
    assert "--num_diffusion_samples=2" in result["command"]
    assert not (bundle / "predictions").exists()


def test_preparation_rejects_unmatched_cohort_fixed_mutations_and_wrong_msa(
    candidate_fixture, tmp_path
):
    path, protocol, _, candidates = candidate_fixture
    path.write_text("".join(json.dumps(row) + "\n" for row in candidates[1:]))
    with pytest.raises(ValueError, match="identical"):
        prepare_jobs(path, protocol, tmp_path / "missing")
    candidates[0]["binder_sequence"] = "W" + "A" * 7
    path.write_text("".join(json.dumps(row) + "\n" for row in candidates))
    with pytest.raises(ValueError, match="fixed binder"):
        prepare_jobs(path, protocol, tmp_path / "mutation")
    candidates[0]["binder_sequence"] = "A" * 8
    path.write_text("".join(json.dumps(row) + "\n" for row in candidates))
    protocol["target_msa"] = "provided"
    msa = tmp_path / "wrong.a3m"
    msa.write_text(">wrong\nAAAA\n")
    with pytest.raises(ValueError, match="query"):
        prepare_jobs(
            path, protocol, tmp_path / "wrong", {hashlib.sha256(b"G" * 12).hexdigest(): str(msa)}
        )


def test_fold_and_pose_are_distinct_and_rigid_invariant(candidate_fixture):
    _, _, reference, _ = candidate_fixture
    target, binder = reference["target_backbone"][:, 1], reference["binder_backbone"][:, 1]
    rotation = np.array([[0, 1, 0], [-1, 0, 0], [0, 0, 1]])
    metrics = geometry_metrics(reference, target @ rotation + 10, binder @ rotation + 10)
    assert metrics["binder_fold_rmsd"] < 1e-10 and metrics["binder_pose_rmsd"] < 1e-10
    displaced = geometry_metrics(reference, target, binder + [0, 0, 6])
    assert displaced["binder_fold_rmsd"] < 1e-10
    assert displaced["binder_pose_rmsd"] == pytest.approx(6)
    with pytest.raises(ValueError, match="Collinear"):
        rigid_fit(
            np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0]]), np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0]])
        )


def test_structural_pipeline_ranks_by_af3_not_desired_pose(candidate_fixture, tmp_path):
    path, protocol, reference, _ = candidate_fixture
    bundle = tmp_path / "bundle"
    plan = prepare_jobs(path, protocol, bundle)
    predictions = fake_outputs(bundle, plan, reference, bad_top_pose=True)
    report = score_bundle(bundle / "manifest.json", predictions, tmp_path / "scores")
    assert report["complete_cohort"] and not report["execution_verified"]
    assert report["arms"]["combined"]["family_pass_rate"]["mean"] == 0
    assert report["arms"]["combined"]["noninferiority_supported"] is None
    rows = [
        json.loads(line) for line in (tmp_path / "scores/candidates.jsonl").read_text().splitlines()
    ]
    assert all(row["selected_sample"].endswith("sample-0") for row in rows)
    assert all(row["metrics"]["pair_iptm"] == 0.9 for row in rows)
    assert all(row["biological_evaluation"] == "not_assessed" for row in rows)


def test_af3_selection_uses_unrounded_ranking_scores(candidate_fixture, tmp_path):
    path, protocol, reference, _ = candidate_fixture
    bundle = tmp_path / "bundle"
    plan = prepare_jobs(path, protocol, bundle)
    predictions = fake_outputs(bundle, plan, reference, bad_top_pose=True)
    for job in plan["jobs"]:
        directory = predictions / job["name"]
        (directory / f"{job['name']}_ranking_scores.csv").write_text(
            "seed,sample,ranking_score\n101,0,0.901\n101,1,0.904\n"
        )
        for sample in range(2):
            summary_path = (
                directory / f"seed-101_sample-{sample}" / "fixture_summary_confidences.json"
            )
            summary = json.loads(summary_path.read_text())
            summary["ranking_score"] = 0.9
            summary_path.write_text(json.dumps(summary))
    report = score_bundle(bundle / "manifest.json", predictions, tmp_path / "scores")
    assert report["complete_cohort"]
    assert report["arms"]["combined"]["family_pass_rate"]["mean"] == 1.0


@pytest.mark.parametrize(
    "fault", ["missing", "nan", "input_changed", "wrong_sequence", "extra_sample"]
)
def test_invalid_outputs_do_not_disappear_from_comparisons(candidate_fixture, tmp_path, fault):
    path, protocol, reference, _ = candidate_fixture
    bundle = tmp_path / "bundle"
    plan = prepare_jobs(path, protocol, bundle)
    predictions = fake_outputs(bundle, plan, reference)
    job = plan["jobs"][0]
    directory = predictions / job["name"]
    first = directory / "seed-101_sample-0"
    if fault == "missing":
        (first / "fixture_model.cif").unlink()
    elif fault == "nan":
        summary = json.loads((first / "fixture_summary_confidences.json").read_text())
        summary["ranking_score"] = float("nan")
        (first / "fixture_summary_confidences.json").write_text(json.dumps(summary))
    elif fault == "input_changed":
        data = json.loads((directory / f"{job['name']}_data.json").read_text())
        data["sequences"][0]["protein"]["templates"] = [{}]
        (directory / f"{job['name']}_data.json").write_text(json.dumps(data))
    elif fault == "extra_sample":
        (directory / "seed-101_sample-2").mkdir()
    else:
        write_prediction(
            first / "fixture_model.cif",
            "G" * 12,
            "W" * 8,
            reference["target_backbone"][:, 1],
            reference["binder_backbone"][:, 1],
        )
    report = score_bundle(bundle / "manifest.json", predictions, tmp_path / "scores")
    assert not report["complete_cohort"]
    assert report["arms"]["frozen_esm3"]["planned_candidates"] == 2
    assert report["arms"]["frozen_esm3"]["valid_output_rate"] == 0.5
    assert "paired_pass_rate_delta_vs_proteinmpnn" not in report["arms"]["combined"]


def test_noninferiority_requires_locked_protocol_and_verified_independent_test_families(
    candidate_fixture,
):
    _, protocol, _, candidates = candidate_fixture
    rows = []
    for i in range(10):
        for candidate in candidates:
            rows.append(
                {
                    **candidate,
                    "candidate_id": f"{i}:{candidate['candidate_id']}",
                    "complex_id": f"complex{i}",
                    "component": f"family{i}",
                    "prediction_status": "complete",
                    "structural_pass": candidate["arm"] != "lora",
                }
            )
    assert (
        summarize_structures(rows, protocol, True)["arms"]["combined"]["noninferiority_supported"]
        is None
    )
    protocol["protocol_state"] = "locked"
    assert (
        summarize_structures(rows, protocol, False)["arms"]["combined"]["noninferiority_supported"]
        is None
    )
    summary = summarize_structures(rows, protocol, True)
    assert summary["arms"]["combined"]["noninferiority_supported"] is True
    assert summary["arms"]["lora"]["noninferiority_supported"] is False
    assert summary["arms"]["combined"]["paired_pass_rate_delta_vs_lora"]["mean"] == 1.0
    for row in rows:
        row["split"] = "application_holdout"
    assert (
        summarize_structures(rows, protocol, True)["arms"]["combined"]["noninferiority_supported"]
        is None
    )


def test_generic_candidate_export_connects_sequence_results_to_af3(candidate_fixture, tmp_path):
    from mertk.candidates import export_candidates
    from mertk.comparison import sequence_metrics

    _, protocol, reference, candidates = candidate_fixture
    first = candidates[0]
    row = {
        "id": first["complex_id"],
        "component": first["component"],
        "split": "test",
        "path": first["reference_npz"],
    }
    results = []
    for arm in protocol["required_arms"]:
        sequences = ["A" * 8, "AG" + "A" * 6]
        results.append(
            {
                "arm": arm,
                "seed": 42,
                "split": "test",
                "smoke_test": False,
                "dataset_identity": first["dataset_identity"],
                "reference_kind": "native",
                "model_identity": {"kind": "fixture-only"},
                "sampling": {"samples": 2, "steps": 2, "temperature": 1.0},
                "cases": [
                    {
                        "id": row["id"],
                        "component": row["component"],
                        "source_sha256": first["reference_sha256"],
                        "mask_sha256": first["mask_sha256"],
                        "sequences": sequences,
                        **sequence_metrics("A" * 8, sequences, reference["design_mask"]),
                    }
                ],
            }
        )
    path = tmp_path / "exported.jsonl"
    exported = export_candidates([row], results, path)
    assert len(exported) == 9
    assert len(prepare_jobs(path, protocol, tmp_path / "exported_bundle")["jobs"]) == 9
    results[0]["smoke_test"] = True
    with pytest.raises(ValueError):
        export_candidates([row], results, tmp_path / "bad.jsonl")
