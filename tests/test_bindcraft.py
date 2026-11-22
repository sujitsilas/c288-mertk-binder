import json

import pytest

from mertk.af3 import prepare_jobs
from mertk.bindcraft import validate_hotspots
from mertk.candidates import import_complexes
from mertk.structures import read_design_chain


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


def import_one(path, output, arm="bindcraft", **kwargs):
    run = output.parent / "run.json"
    run.write_text(json.dumps({"tool_revision": "fixture-only"}))
    return import_complexes(
        [path],
        output,
        arm=arm,
        campaign_id="pilot",
        epitope_id="patch1",
        target_reference=path,
        run_record=run,
        seed=42,
        **kwargs,
    )


def test_import_to_af3_keeps_sequences_and_provenance(complex_pdb, tmp_path):
    rows = import_one(complex_pdb, tmp_path / "bindcraft")
    rows += import_one(complex_pdb, tmp_path / "rfd4", arm="rfd4_proteina")
    assert rows[0]["binder_sequence"] == "A" * 16
    assert len(rows[0]["binder_mapping"]) == 16
    candidates = tmp_path / "candidates.jsonl"
    candidates.write_text("".join(json.dumps(r) + "\n" for r in rows))
    from pathlib import Path

    protocol = json.loads(Path("configs/af3_evaluation.json").read_text())
    protocol["target_msa"] = "single_sequence"
    plan = prepare_jobs(candidates, protocol, tmp_path / "bundle")
    assert len(plan["jobs"]) == 2
    assert all(j["candidate"]["biological_evaluation"] == "not_assessed" for j in plan["jobs"])


def test_hotspots_and_numbering_are_validated(complex_pdb):
    _, mapping = read_design_chain(complex_pdb, "A")
    assert validate_hotspots("A1,A3-5", "A", mapping) == "A1,A3-5"
    with pytest.raises(ValueError, match="absent"):
        validate_hotspots("A999", "A", mapping)
    with pytest.raises(ValueError, match="selected chain"):
        validate_hotspots("B5", "A", mapping)


def test_import_rejects_invalid_chains_and_existing_output(complex_pdb, tmp_path):
    with pytest.raises(ValueError, match="different chains"):
        import_one(complex_pdb, tmp_path / "bad", binder_chain="A")
    with pytest.raises(ValueError, match="absent"):
        read_design_chain(complex_pdb, "Z")
    import_one(complex_pdb, tmp_path / "ok")
    with pytest.raises(FileExistsError):
        import_one(complex_pdb, tmp_path / "ok")
