"""Synthetic safety cases for gating raw generators, not experimental binder tests."""

import runpy

import numpy as np
import pytest
from Bio.PDB import PDBIO, Atom, Chain, Model, Residue, Structure


def write(path, target, binder=None, glycan=None):
    model = Model.Model(0)
    structure = Structure.Structure("synthetic")
    structure.add(model)
    for chain_id, residue_name, points in [("A", "GLY", target), ("B", "ALA", binder)]:
        if points is None:
            continue
        chain = Chain.Chain(chain_id)
        model.add(chain)
        for i, point in enumerate(points, 1):
            residue = Residue.Residue((" ", i, " "), residue_name, " ")
            chain.add(residue)
            for name, shift in [
                ("N", [-1, 1, 0]),
                ("CA", [0, 0, 0]),
                ("C", [1, 0, 0]),
                ("O", [1, 1, 0]),
            ]:
                residue.add(Atom.Atom(name, point + shift, 0, 1, " ", name, i, element=name[0]))
    if glycan is not None:
        chain = Chain.Chain("G")
        model.add(chain)
        residue = Residue.Residue(("H_NAG", 1, " "), "NAG", " ")
        chain.add(residue)
        residue.add(Atom.Atom("C1", glycan, 0, 1, " ", "C1", 1, element="C"))
    writer = PDBIO()
    writer.set_structure(structure)
    writer.save(str(path))


@pytest.fixture
def inputs(tmp_path):
    target = np.array([[0.0, 0, 0], [3.8, 0, 0], [3.8, 3.8, 0], [3.8, 3.8, 3.8]])
    binder = target[:3] + [0, 0, 5]
    reference, output = tmp_path / "reference.pdb", tmp_path / "design.pdb"
    write(reference, target)
    write(output, target, binder)
    analyze = runpy.run_path("scripts/hpc/check_course_outputs.py")["analyze"]
    return target, binder, reference, output, analyze


def test_geometry_pass_and_exact_chain_assignment(inputs):
    _, _, reference, output, analyze = inputs
    row = analyze(output, reference)
    assert row["preliminary_geometry_pass"]
    assert row["target_chain"] == "A" and row["binder_chain"] == "B"
    assert row["target_ca_rmsd_A"] == pytest.approx(0, abs=1e-7)


def test_target_distortion_cannot_be_hidden_by_context_restoration(inputs):
    target, binder, reference, output, analyze = inputs
    broken = target.copy()
    broken[-1] += [30, 0, 0]
    write(output, broken, binder)
    before = output.read_bytes()
    row = analyze(output, reference)
    assert not row["preliminary_geometry_pass"]
    assert any("RMSD" in message for message in row["rejection_reasons"])
    assert output.read_bytes() == before


def test_collapsed_binder_backbone_rejected(inputs):
    target, binder, reference, output, analyze = inputs
    binder[1] = binder[0] + [0.4, 0, 0]
    write(output, target, binder)
    row = analyze(output, reference)
    assert not row["preliminary_geometry_pass"]
    assert row["adjacent_ca"]["binder"]["outside_2p7_to_4p5_A"] > 0


def test_observed_glycan_clash_rejected_in_native_context(inputs, tmp_path):
    target, binder, reference, output, analyze = inputs
    context = tmp_path / "context.pdb"
    write(context, target, glycan=binder[0])
    row = analyze(output, reference, context_path=context, hotspots=(1, 2, 3))
    assert not row["preliminary_geometry_pass"]
    assert row["observed_glycan_clearance_A"] == pytest.approx(0, abs=1e-7)
