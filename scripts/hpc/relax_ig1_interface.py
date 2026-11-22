"""Test whether fixed-pose binder side-chain relaxation resolves native Ig1 clashes."""

import argparse
import json
from pathlib import Path

import numpy as np
from check_course_outputs import analyze, structure

from mertk.structural import rigid_fit

BACKBONE = {"N", "CA", "C", "O", "OXT"}


def placed_complex(start, receptor, output):
    original = structure(start)
    native = structure(receptor)["A"]
    receptor_ca = np.array([r["CA"].coord for r in native])
    predicted_ca = np.array([r["CA"].coord for r in original["A"]])
    rotation, shift = rigid_fit(predicted_ca, receptor_ca)
    for a in original["B"].get_atoms():
        a.coord = a.coord @ rotation + shift
    from Bio.PDB import PDBIO

    original.detach_child("A")
    native = native.copy()
    original.add(native)
    # Preserve receptor A then binder B as required by the interface analysis.
    ligand = original["B"]
    original.detach_child("B")
    original.add(ligand)
    writer = PDBIO()
    writer.set_structure(original)
    writer.save(str(output))


def closest_pair(path, backbone_only=False):
    model = structure(path)
    a, b = (
        [
            atom
            for atom in model[chain].get_atoms()
            if atom.element not in {"H", "D"} and (not backbone_only or atom.name in BACKBONE)
        ]
        for chain in ("A", "B")
    )
    x, y = (np.array([atom.coord for atom in atoms]) for atoms in (a, b))
    distances = np.linalg.norm(x[:, None] - y[None, :], axis=-1)
    i, j = np.unravel_index(distances.argmin(), distances.shape)
    return {
        "distance_A": float(distances[i, j]),
        "atoms": [f"{atom.parent.resname}{atom.parent.id[1]}:{atom.name}" for atom in (a[i], b[j])],
        "pairs_below_2A": int((distances < 2).sum()),
    }


def main(suite):
    import pyrosetta as pr

    pr.init("-mute all -constant_seed -jran 112026")
    manifest = json.loads((suite / "manifest.json").read_text())
    reference = suite / "target.pdb"
    context = suite / "37NT.cif"
    results = {}
    for key, hotspots in manifest.items():
        print("Relaxing", key, flush=True)
        before = suite / f"{key}_placed.pdb"
        after = suite / f"{key}_relaxed.pdb"
        placed_complex(suite / f"{key}.pdb", reference, before)
        pose = pr.pose_from_pdb(str(before))
        move = pr.rosetta.core.kinematics.MoveMap()
        move.set_bb(False)
        move.set_jump(False)
        move.set_chi(False)
        for index in range(1, pose.total_residue() + 1):
            move.set_chi(index, pose.pdb_info().chain(index) == "B")
        relax = pr.rosetta.protocols.relax.FastRelax(pr.get_fa_scorefxn(), 1)
        relax.set_movemap(move)
        relax.set_movemap_disables_packing_of_fixed_chi_positions(True)
        relax.max_iter(100)
        relax.apply(pose)
        pose.dump_pdb(str(after))
        for label, path in (("before", before), ("after", after)):
            report = analyze(path, reference, context_path=context, hotspots=hotspots)
            report["closest_pair"] = closest_pair(path)
            report["closest_backbone_pair"] = closest_pair(path, backbone_only=True)
            results.setdefault(key, {})[label] = report
        (suite / "assessment.json").write_text(json.dumps(results, indent=2) + "\n")
        print(
            key,
            results[key]["before"]["closest_pair"],
            results[key]["after"]["closest_pair"],
            flush=True,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("suite", type=Path)
    main(parser.parse_args().suite)
