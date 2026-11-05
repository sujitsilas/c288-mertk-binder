import json
from pathlib import Path

import numpy as np
import pytest
import torch

from mertk.atom_decoder import AtomDecoder
from mertk.bindcraft import (
    export_candidates,
    import_complexes,
    prepare_job,
    read_design_chain,
    validate_hotspots,
)
from mertk.dataset import read_records
from mertk.training import audit_splits, file_hash


@pytest.fixture
def complex_pdb(tmp_path):
    path = tmp_path / "complex.pdb"
    lines, serial = [], 1
    for chain, y in [("A", 0), ("B", 5)]:
        for i in range(1, 17):
            for atom, dx, dy, element in [
                ("N", -1, 1, "N"),
                ("CA", 0, 0, "C"),
                ("C", 1, 0, "C"),
                ("O", 1, 1, "O"),
                ("CB", 0, -1, "C"),
            ]:
                lines.append(
                    f"ATOM  {serial:5d} {atom:>4s} ALA {chain}{i:4d}    {i * 3.8 + dx:8.3f}{y + dy:8.3f}{0:8.3f}  1.00 80.00          {element:>2}\n"
                )
                serial += 1
        lines.append("TER\n")
    path.write_text("".join(lines) + "END\n")
    return path


def test_import_preserves_original_sequence_and_marks_application_holdout(complex_pdb, tmp_path):
    output = tmp_path / "imported"
    rows = import_complexes([complex_pdb], output)
    assert rows[0]["split"] == "application_holdout"
    assert "AAAAAAAAAAAAAAAA" in (output / "bindcraft_original.fasta").read_text()
    with np.load(output / rows[0]["path"]) as example:
        assert example["design_mask"].all()
        assert example["binder_oxygen"].shape == (16, 3)
    provenance = json.loads((output / "provenance.json").read_text())
    assert provenance[0]["binder_chain"] == "B"
    assert len(provenance[0]["binder_mapping"]) == 16


def test_hotspots_and_length_budget_are_validated(complex_pdb, tmp_path):
    _, mapping = read_design_chain(complex_pdb, "A")
    assert validate_hotspots("A1,A3-5", "A", mapping) == "A1,A3-5"
    with pytest.raises(ValueError, match="absent"):
        validate_hotspots("A999", "A", mapping)
    with pytest.raises(ValueError, match="selected chain"):
        validate_hotspots("B5", "A", mapping)
    with pytest.raises(ValueError, match="token budget"):
        prepare_job(complex_pdb, "A", tmp_path / "job", lengths=(500, 501))


def test_import_rejects_same_chain_or_missing_chain(complex_pdb, tmp_path):
    with pytest.raises(ValueError, match="different chains"):
        import_complexes([complex_pdb], tmp_path / "bad", binder_chain="A")
    with pytest.raises(ValueError, match="absent"):
        read_design_chain(complex_pdb, "Z")


def test_candidate_export_keeps_original_control_and_all_redesigns(complex_pdb, tmp_path):
    output = tmp_path / "imported"
    import_complexes([complex_pdb], output)
    rows = read_records(output / "examples.jsonl")
    results = [
        {
            "reference_kind": "bindcraft_seed",
            "arm": "combined",
            "seed": 42,
            "cases": [{"id": rows[0]["id"], "sequences": ["G" + "A" * 15, "A" * 15 + "G"]}],
        }
    ]
    candidates = export_candidates(rows, results, tmp_path / "candidates.jsonl")
    assert len(candidates) == 3
    assert candidates[0]["arm"] == "bindcraft_original"
    assert candidates[0]["binder_sequence"] == "A" * 16
    assert all(c["structure_evaluation"] == "pending" for c in candidates)
    complex_pdb.write_text(complex_pdb.read_text() + "REMARK modified\n")
    with pytest.raises(ValueError, match="changed after import"):
        export_candidates(rows, results, tmp_path / "bad.jsonl")


def test_application_sequences_use_same_split_audit_keys(complex_pdb, tmp_path):
    rows = import_complexes([complex_pdb], tmp_path / "imported")
    assert rows[0]["target_key"].startswith("seq_")
    training_row = dict(rows[0], split="train", component="train_family", pdb_id="TRAIN")
    with pytest.raises(ValueError, match="leakage"):
        audit_splits([training_row, rows[0]])


def test_application_inference_keeps_training_manifest_and_uses_seed_retention(
    complex_pdb, tmp_path
):
    from mertk.comparison import run_comparison

    output = tmp_path / "imported"
    import_complexes([complex_pdb], output)
    rows = read_records(output / "examples.jsonl")
    path = Path(rows[0]["path"])
    with np.load(path) as data:
        example = dict(data)
    for chain in ("binder", "target"):
        example[f"{chain}_structure_tokens"] = np.full(16, 4096, dtype=np.int64)
    np.savez_compressed(path, **example)
    training_manifest = tmp_path / "training.jsonl"
    # A checkpoint fixture references a disjoint training manifest, never this holdout.
    training_manifest.write_text(
        json.dumps(
            {
                "id": "train",
                "split": "train",
                "component": "train",
                "pdb_id": "TRAIN",
                "binder_key": "train_binder",
                "target_key": "train_target",
            }
        )
        + "\n"
    )
    config = {
        "model_type": "atom_decoder",
        "hidden": 16,
        "manifest": str(training_manifest),
        "device": "cpu",
        "seed": 42,
        "precision": "fp32",
    }
    checkpoint = tmp_path / "fixture.pt"
    torch.save(
        {
            "metadata": {
                "config": config,
                "smoke_test": False,
                "weights": {"kind": "atom_decoder"},
                "manifest_sha256": file_hash(training_manifest),
            },
            "adapter": AtomDecoder(16).adapter_state_dict(),
        },
        checkpoint,
    )
    result = run_comparison(
        config,
        "atom_decoder",
        tmp_path / "result.json",
        checkpoint,
        samples=2,
        sampling_steps=2,
        split="application_holdout",
        application_manifest=output / "examples.jsonl",
    )
    assert result["reference_kind"] == "bindcraft_seed"
    assert "mean_seed_retention" in result["cases"][0]
    assert "mean_recovery" not in result["cases"][0]
