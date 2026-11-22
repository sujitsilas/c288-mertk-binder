import json
import runpy
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pytest
from Bio.PDB import MMCIFIO, PDBIO, Atom, Chain, Model, PDBParser, Residue, Structure

pytest.importorskip("matplotlib")


def write_structure(path, chains):
    structure = Structure.Structure("synthetic_plot_test")
    model = Model.Model(0)
    structure.add(model)
    for key, name, points in chains:
        chain = Chain.Chain(key)
        model.add(chain)
        for i, point in enumerate(points, 1):
            residue = Residue.Residue((" ", i, " "), name, " ")
            for atom, offset in (
                ("N", [-1, 1, 0]),
                ("CA", [0, 0, 0]),
                ("C", [1, 0, 0]),
                ("O", [1, 1, 0]),
            ):
                residue.add(Atom.Atom(atom, point + offset, 0, 1, " ", atom, i, element=atom[0]))
            chain.add(residue)
    writer = PDBIO()
    writer.set_structure(structure)
    writer.save(str(path))


def test_plot_assignment_geometry_and_duplicate_accounting(tmp_path):
    (tmp_path / "inputs").mkdir()
    designs = tmp_path / "bindcraft/designs"
    designs.mkdir(parents=True)
    target = np.array([[0.0, 0, 0], [4, 0, 0], [4, 4, 0], [4, 4, 4]])
    binder = target[:3] + [0, 0, 5]
    write_structure(tmp_path / "inputs/target.pdb", [("A", "GLY", target)])
    # Target deliberately has a different chain ID and translated coordinates.
    write_structure(
        designs / "generated.pdb", [("Z", "GLY", target + 30), ("A", "ALA", binder + 30)]
    )
    (designs / "copy.pdb").write_bytes((designs / "generated.pdb").read_bytes())
    write_structure(designs / "wrong_target.pdb", [("Z", "VAL", target), ("A", "ALA", binder)])
    collect = runpy.run_path("scripts/plot_pilot.py")["collect"]
    rows, skipped, _ = collect(tmp_path)
    assert len(rows) == 1
    assert rows[0]["target_chain"] == "Z"
    assert rows[0]["binder_chain"] == "A"
    assert rows[0]["target_ca_rmsd"] == pytest.approx(0, abs=1e-8)
    expected = np.linalg.norm(target[:, None] - binder[None, :], axis=-1)
    assert rows[0]["ca_contact_pairs"] == int((expected <= 8).sum())
    assert len(skipped) == 2
    subprocess.run([sys.executable, "scripts/plot_pilot.py", str(tmp_path)], check=True)
    assert len(list((tmp_path / "figures").glob("*.png"))) == 3
    assert len(json.loads((tmp_path / "figures/geometry.json").read_text())["records"]) == 1


def test_no_designs_produces_no_result_figures(tmp_path):
    (tmp_path / "inputs").mkdir()
    write_structure(tmp_path / "inputs/target.pdb", [("A", "GLY", np.eye(3))])
    subprocess.run([sys.executable, "scripts/plot_pilot.py", str(tmp_path)], check=True)
    assert not list((tmp_path / "figures").glob("*.png"))
    assert "No complete" in (tmp_path / "figures/README.md").read_text()


def test_assessment_cohort_uses_both_arms_and_one_complex_per_seed(tmp_path, monkeypatch):
    (tmp_path / "inputs").mkdir()
    accepted = tmp_path / "bindcraft/designs/Accepted"
    accepted.mkdir(parents=True)
    generated = tmp_path / "rfd4_proteina/inference_samples"
    generated.mkdir(parents=True)
    target = np.array([[0.0, 0, 0], [4, 0, 0], [4, 4, 0], [4, 4, 4]])
    binder = target[:3] + [0, 0, 5]
    write_structure(tmp_path / "inputs/target.pdb", [("A", "GLY", target)])
    source = accepted / "fixture_l3_s17_mpnn1.pdb"
    write_structure(source, [("Z", "GLY", target), ("A", "ALA", binder)])
    # Distinct coordinates, same seed: only one BindCraft trajectory may be nominated.
    write_structure(
        accepted / "fixture_l3_s17_mpnn2.pdb", [("Z", "GLY", target), ("A", "ALA", binder + 1)]
    )
    writer = MMCIFIO()
    writer.set_structure(PDBParser(QUIET=True).get_structure("fixture", str(source)))
    writer.save(str(generated / "fixture_generated.cif"))
    # Correct sequence and close binder contacts cannot rescue a deformed target.
    deformed = target.copy()
    deformed[-1] += [50, 0, 0]
    distorted = tmp_path / "distorted.pdb"
    write_structure(distorted, [("Z", "GLY", deformed), ("A", "ALA", binder)])
    writer.set_structure(PDBParser(QUIET=True).get_structure("distorted", str(distorted)))
    writer.save(str(generated / "distorted_generated.cif"))
    (tmp_path / "manifest.json").write_text(
        json.dumps(
            {
                "campaign_id": "synthetic_test",
                "rfd4_seed": 101,
                "tools": {
                    arm: {"revision": "fixture-only"} for arm in ("bindcraft", "rfd4_proteina")
                },
                "jobs": {
                    f"{arm}_generate": {"job_id": "0"} for arm in ("bindcraft", "rfd4_proteina")
                },
            }
        )
    )
    main = runpy.run_path("scripts/hpc/prepare_pilot_assessment.py")["main"]
    monkeypatch.setattr(sys, "argv", ["prepare_pilot_assessment.py", str(tmp_path)])
    monkeypatch.setattr(
        subprocess, "run", lambda *a, **k: SimpleNamespace(stdout="fixture", stderr="")
    )
    main()
    selected = json.loads((tmp_path / "assessment/selection.json").read_text())
    assert selected["counts"] == {"bindcraft": 1, "rfd4_proteina": 1}
    assert len(selected["rejected"]) == 1
    assert "Target C-alpha RMSD" in selected["rejected"][0]["reason"]
    plan = json.loads((tmp_path / "assessment/af3/manifest.json").read_text())
    assert len(plan["jobs"]) == 2
    assert plan["protocol"]["target_msa"] == "single_sequence"
    assert plan["protocol"]["protocol_state"] == "locked"
