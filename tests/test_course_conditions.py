"""A parser's residue renumbering must not hide changed conditions or geometry."""

import runpy

import numpy as np
import pytest

audit = runpy.run_path("scripts/hpc/course_condition_checks.py")["audit_arrays"]


@pytest.fixture
def arrays():
    names = ["N", "CA", "C", "O", "CB", "N", "CA", "C", "O", "CB", "N"]
    points = np.array(
        [
            [0, 0, 0],
            [1, 1, 0],
            [2, 0, 0],
            [2, -1, 0],
            [1, 1, 1],
            [3, 0, 0],
            [4, 1, 0],
            [5, 0, 1],
            [5, -1, 1],
            [4, 1, 1],
            [0, 0, 5],
        ],
        dtype=float,
    )
    target = np.array([True] * 10 + [False])
    backbone = target & np.isin(names, ["N", "CA", "C", "O"])
    return dict(
        chain_id=np.array(["A"] * 10 + ["B"]),
        res_id=np.array([1] * 5 + [2] * 5 + [1]),
        res_name=np.array(["ALA"] * 5 + ["SER"] * 5 + ["GLY"]),
        atom_name=np.array(names),
        coord=points,
        known=backbone,
        sequence=target,
        hot=np.array([False] * 4 + [True] + [False] * 6),
        reference={
            "residues": [
                {
                    "author_number": n,
                    "residue_name": name,
                    "backbone": {
                        atom: points[offset + i].tolist()
                        for i, atom in enumerate(["N", "CA", "C", "O"])
                    },
                }
                for n, name, offset in [(178, "ALA", 0), (183, "SER", 5)]
            ],
            "hotspot_atoms": [[178, "CB"]],
        },
    )


def test_renumbered_native_crop_passes_with_explicit_map(arrays):
    row = audit(**arrays)
    assert row["parsed_to_author_residue_map"] == {1: 178, 2: 183}
    assert row["hotspot_residues"] == [178]


def test_rigid_transformation_is_allowed(arrays):
    arrays["coord"] = arrays["coord"] @ np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]]) + 20
    assert audit(**arrays)["max_backbone_deviation_A"] < 1e-6


def test_wrong_hotspot_with_same_atom_count_rejected(arrays):
    arrays["hot"][[4, 9]] = [False, True]
    with pytest.raises(ValueError, match="Hotspot atom identities"):
        audit(**arrays)


def test_changed_sequence_order_rejected(arrays):
    arrays["res_name"][:5] = "SER"
    arrays["res_name"][5:10] = "ALA"
    with pytest.raises(ValueError, match="identities/order"):
        audit(**arrays)


def test_distorted_native_coordinates_rejected(arrays):
    arrays["coord"][1] += 5
    with pytest.raises(ValueError, match="geometry changed"):
        audit(**arrays)


def test_wrong_coordinate_condition_mask_rejected(arrays):
    arrays["known"][4] = True
    with pytest.raises(ValueError, match="coordinate mask"):
        audit(**arrays)
