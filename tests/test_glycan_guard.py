import numpy as np
import pytest

from mertk.glycan_guard import frames, losses, prepare_constraints


def example():
    xyz = np.zeros((3, 37, 3))
    xyz[0, :3] = [[0, 1, 0], [0, 0, 0], [1, 0, 0]]
    xyz[1, :3] = [[40, 1, 0], [40, 0, 0], [41, 0, 0]]
    xyz[2, :] = [0.4, 0, 3]
    mask = np.zeros((3, 37))
    mask[:, :3] = 1
    spec = {
        "target_length": 2,
        "anchor_indices": [0],
        "sequon_indices": [0],
        "glycan_local_xyz": [[0, 0, 3]],
        "glycan_atom_sites": [0],
        "heavy_atom_clearance_A": 3.5,
        "sequon_ca_exclusion_A": 15.0,
    }
    return xyz, mask, spec


def test_collision_and_remote_binder():
    xyz, mask, spec = example()
    assert min(losses(xyz, mask, spec).values()) > 0.5
    xyz[2] += 100
    assert losses(xyz, mask, spec) == {"glycan_volume": 0, "sequon_volume": 0}


def test_rigid_body_invariance():
    xyz, mask, spec = example()
    q, _ = np.linalg.qr(np.random.default_rng(52).normal(size=(3, 3)))
    assert np.linalg.det(q) > 0
    moved = xyz @ q + np.array([14, -27, 3])
    assert losses(moved, mask, spec) == pytest.approx(losses(xyz, mask, spec), abs=1e-6)
    f = frames(xyz[0, :3])
    np.testing.assert_allclose(f @ f.T, np.eye(3), atol=1e-7)


def test_nonexistent_atoms_do_not_clash():
    xyz, mask, spec = example()
    xyz[2, :3] += 100
    assert losses(xyz, mask, spec) == {"glycan_volume": 0, "sequon_volume": 0}


def test_native_attachments():
    from pathlib import Path

    crop = Path("outputs/design/ig_only_20261113/ig_pair.pdb")
    cif = Path("data/references/37NT/37NT.cif")
    if not crop.exists() or not cif.exists():
        pytest.skip("Native experimental inputs not installed")
    spec = prepare_constraints(crop, cif)
    assert [a["native_asn"] for a in spec["attachments"]] == [114, 170, 207, 215, 234]
    assert spec["target_length"] == 193
    assert all(a["covalent_distance_A"] < 1.6 for a in spec["attachments"])
