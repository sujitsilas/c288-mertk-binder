"""Viewer transforms must preserve the scientific meaning of native overlays."""

import io
import json

import numpy as np
import pytest
from Bio.PDB import Atom, Chain, PDBParser, Residue
from Bio.SeqUtils import seq3

from mertk.prediction_review import native_numbering, nearest_points, pdb_chain, safe_json


def chain_fixture():
    chain = Chain.Chain("A")
    sequence = "ACDEFGHIKLMNPQRSTVWYACDEFHIKLMNPQ"
    for index, aa in enumerate(sequence, 1):
        residue = Residue.Residue((" ", index, " "), seq3(aa).upper(), " ")
        for j, (name, element) in enumerate([("N", "N"), ("CA", "C"), ("C", "C"), ("H", "H")]):
            residue.add(
                Atom.Atom(
                    name,
                    np.array([index, j, index % 3], dtype=float),
                    80.0,
                    1.0,
                    " ",
                    name,
                    index * 4 + j,
                    element=element,
                )
            )
        chain.add(residue)
    return chain, sequence


def test_generator_numbering_is_mapped_by_exact_sequence():
    chain, seq = chain_fixture()
    mapping = native_numbering(chain, "G" * 92 + seq + "V" * 10)
    assert mapping[(" ", 1, " ")] == 93
    assert mapping[(" ", len(seq), " ")] == 92 + len(seq)


def test_ambiguous_sequence_does_not_get_an_arbitrary_native_offset():
    chain, seq = chain_fixture()
    with pytest.raises(ValueError, match="unique"):
        native_numbering(chain, seq + "G" * 10 + seq)


def test_rigid_overlay_and_renumbering_do_not_mutate_prediction():
    chain, seq = chain_fixture()
    before = np.array([a.coord.copy() for a in chain.get_atoms()])
    mapping = native_numbering(chain, "G" * 92 + seq)
    rot = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    shift = np.array([4, 8, -3])
    text = pdb_chain(chain, mapping, (rot, shift))
    displayed = PDBParser(QUIET=True).get_structure("view", io.StringIO(text))[0]["A"]
    assert list(displayed)[0].id[1] == 93
    assert np.allclose(displayed[93]["CA"].coord, chain[1]["CA"].coord @ rot + shift)
    assert np.array_equal(before, np.array([a.coord for a in chain.get_atoms()]))
    assert list(chain)[0].id[1] == 1
    assert not any(a.element == "H" for a in displayed.get_atoms())


def test_chunked_clash_distance_matches_direct_calculation():
    rng = np.random.default_rng(8)
    x, y = rng.normal(size=(260, 3)), rng.normal(size=(42, 3))
    distance, index = nearest_points(x, y)
    direct = np.linalg.norm(x[:, None] - y[None], axis=-1)
    assert np.allclose(distance, direct.min(axis=1))
    assert np.array_equal(index, direct.argmin(axis=1))


def test_embedded_data_cannot_close_the_json_script_element():
    value = {"name": "</script><script>alert(1)</script>"}
    encoded = safe_json(value)
    assert "<" not in encoded
    assert json.loads(encoded) == value
