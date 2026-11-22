import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest
from Bio.PDB import Chain, PDBParser, Residue

from mertk.cleavage_protector import _map_target, rotate_psi
from mertk.glycan_guard import losses
from mertk.protector_screen import disable_score_cutoffs
from mertk.provenance import file_hash
from mertk.target_policy import load_policy, require_active_target


def module(name):
    scripts = Path("scripts/hpc").resolve()
    sys.path.insert(0, str(scripts))
    try:
        spec = importlib.util.spec_from_file_location(name, scripts / f"{name}.py")
        result = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(result)
        return result
    finally:
        sys.path.pop(0)


def test_protector_scope_does_not_reopen_retired_patch():
    policy = load_policy()
    record = {"workflow": "cleavage_protector_v1", "aim": "aim2_shedding_protector"}
    assert require_active_target("Stalk_P490S491", [488, 490, 491], record, policy) == "Stalk"
    for name, hotspots in [("FNIIIb_F472", [472]), ("Stalk_P490S491", [472, 490, 491])]:
        with pytest.raises(ValueError):
            require_active_target(name, hotspots, record, policy)
    with pytest.raises(ValueError):
        require_active_target("Stalk_P490S491", [488, 490, 491], {}, policy)


def test_stalk_mapping_ignores_interleaved_hetero_residues_and_rejects_missing_residues():
    reference = Chain.Chain("A")
    target = Chain.Chain("X")
    for index in range(1, 142):
        name = "PRO" if index == 105 else "ALA"
        reference.add(Residue.Residue((" ", 385 + index, " "), name, ""))
        target.add(Residue.Residue((" ", index, " "), name, ""))
        if index == 50:
            target.add(Residue.Residue(("W", index, " "), "HOH", ""))
    mapped = _map_target(target, reference)
    assert list(mapped) == list(range(386, 527))
    assert mapped[490].resname == "PRO"
    assert mapped[491] is target[106]

    target.detach_child((" ", 50, " "))
    with pytest.raises(ValueError, match="complete protector target sequence"):
        _map_target(target, reference)
    reference.detach_child((" ", 435, " "))
    with pytest.raises(ValueError, match="complete native"):
        _map_target(target, reference)


def test_membrane_loss_is_masked_and_rigid_body_invariant():
    p = np.zeros((4, 37, 3))
    for i in range(3):
        p[i, :3] = np.array([[0, 1, 0], [0, 0, 0], [1, 0, 0]]) + [0, 0, i * 5]
    p[3, :3] = [30, 0, 1]
    m = np.zeros((4, 37))
    m[:, :3] = 1
    s = dict(
        target_length=3,
        anchor_indices=[0],
        sequon_indices=[0],
        glycan_local_xyz=[[0, 0, 3]],
        glycan_atom_sites=[0],
        heavy_atom_clearance_A=3.5,
        sequon_ca_exclusion_A=15,
        membrane={
            "outer_indices": [0],
            "inner_indices": [2],
            "origin_index": 0,
            "extracellular_margin_A": 5,
        },
    )
    assert losses(p, m, s)["membrane_volume"] > 0
    r, _ = np.linalg.qr(np.random.default_rng(4).normal(size=(3, 3)))
    assert losses(p @ r + [3, 7, 2], m, s) == pytest.approx(losses(p, m, s), abs=1e-6)
    p[3, :3, 2] = -20
    assert losses(p, m, s)["membrane_volume"] == 0


def test_torsion_perturbation_preserves_covalent_geometry():
    source = Path("outputs/design/cleavage_protector_20261117/AF_Q12866.pdb")
    if not source.exists():
        pytest.skip("AFDB model not installed")
    model = PDBParser(QUIET=True).get_structure("af", source)[0]
    rotated = rotate_psi(model, 490, 25)
    assert rotated["A"][490]["C"] - rotated["A"][491]["N"] == pytest.approx(
        model["A"][490]["C"] - model["A"][491]["N"], abs=1e-5
    )
    assert np.linalg.norm(rotated["A"][491]["CA"].coord - model["A"][491]["CA"].coord) > 0.1
    np.testing.assert_array_equal(rotated["A"][480]["CA"].coord, model["A"][480]["CA"].coord)


def test_budget_and_held_assessment():
    submit = module("submit_protector_round")
    record = json.loads(Path("configs/cleavage_protector_round.json").read_text())
    rows = submit.job_specs(record)
    assert sum(submit.seconds(r[6]) for r in rows if r[5]) == 8400
    assessment = next(r for r in rows if r[0] == "af3_assessment")
    assert assessment[7] is True
    assert assessment[3] == ["protector_select"]


def test_empty_protector_screen_skips_assessment(tmp_path):
    run = module("protector_round")
    (tmp_path / "manifest.json").write_text(
        json.dumps(
            {"experiments": {"a": {}}, "af3_protocol": {}, "canonical_target_sequence": "AAA"}
        )
    )
    run.select_parents(tmp_path, manage_scheduler=False)
    assert json.loads((tmp_path / "assessment/manifest.json").read_text())["cases"] == []
    assert (tmp_path / "assessment/execution_skipped.json").exists()


def test_exploratory_filter_keeps_scalar_and_residue_metrics():
    settings = {
        "1_i_pAE": {"threshold": 0.35, "higher": False},
        "Average_InterfaceAAs": {"Y": {"threshold": 1, "higher": True}},
    }
    result = disable_score_cutoffs(settings)
    assert result["1_i_pAE"] == {"threshold": None, "higher": False}
    assert result["Average_InterfaceAAs"]["Y"] == {"threshold": None, "higher": True}


@pytest.mark.parametrize("mode, count", [("exploratory", 1), ("original", 0)])
def test_ranked_stalk_nomination_without_strict_pass(tmp_path, mode, count):
    run = module("protector_round")
    record = dict(
        experiments={"pilot": {}},
        af3_protocol={},
        canonical_target_sequence="AAA",
        protector_screen=mode,
    )
    (tmp_path / "manifest.json").write_text(json.dumps(record))
    campaign = tmp_path / "campaigns/pilot"
    campaign.mkdir(parents=True)
    records = []
    for name, stage, distances in [
        ("early", "Clash flagged", [0.1, 0.2]),
        ("further", "Scored backbone", [3.2, 4.0]),
        ("closer", "Scored backbone", [3.7, 3.8]),
    ]:
        parent = campaign / f"{name}.pdb"
        parent.write_text(name)
        records.append(
            dict(
                name=name,
                stage=stage,
                cleavage_atom_min_distances_A=distances,
                iptm=None,
                coordinate_file=parent.name,
                sha256=file_hash(parent),
                binder_sequence="GG",
                preliminary_model_pass=False,
            )
        )
    (campaign / "protector_quality.json").write_text(
        json.dumps(dict(records=records, accepted_variants=0, preliminary_model_passes=0))
    )
    run.select_parents(tmp_path, manage_scheduler=False)
    cases = json.loads((tmp_path / "assessment/manifest.json").read_text())["cases"]
    assert len(cases) == count
    if cases:
        assert cases[0]["parent"] == "campaigns/pilot/closer.pdb"
