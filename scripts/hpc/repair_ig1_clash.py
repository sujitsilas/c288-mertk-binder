"""Diagnose and repair the MerTK Lys112 contact with the original Ig1 mpnn8."""

import argparse
import json
from pathlib import Path

import numpy as np
from check_course_outputs import analyze, structure
from relax_ig1_interface import closest_pair
from score_ig1_refolds import refold, write_json

from mertk.provenance import file_hash
from mertk.structural import rigid_fit, rmsd
from mertk.structures import read_design_chain

HOTSPOTS = (122, 148)
GEOMETRY_CASES = ("lys112_rotamer", "i32g", "local_backbone", "i32g_local_backbone")
MULTISITE_CASES = {
    "local_backbone_multisite": {},
    "i32g_w83f": {32: "GLY", 83: "PHE"},
    "i32g_w83l": {32: "GLY", 83: "LEU"},
    "i32g_w83f_k15q": {32: "GLY", 83: "PHE", 15: "GLN"},
}


def contacts(path, reference):
    model = structure(path)
    ligand = np.array([a.coord for a in model["B"].get_atoms() if a.element not in {"H", "D"}])
    target = [r for r in model["A"] if r.id[0] == " "]
    positions = [r.id[1] for r in structure(reference)["A"] if r.id[0] == " "]
    if len(positions) != len(target):
        raise ValueError("MerTK residue counts differ")
    return {
        n
        for n, residue in zip(positions, target, strict=True)
        if any(
            np.linalg.norm(atom.coord - ligand, axis=1).min() <= 4.5
            for atom in residue
            if atom.element not in {"H", "D"}
        )
    }


def displacement(path, start):
    target, ligand = [
        [read_design_chain(p, chain)[0].backbone[:, 1] for p in (start, path)]
        for chain in ("A", "B")
    ]
    r, t = rigid_fit(target[1], target[0])
    return rmsd(ligand[1] @ r + t, ligand[0])


def review(path, suite, start, baseline, max_rmsd):
    reference = suite / "target.pdb"
    row = analyze(path, reference, context_path=suite / "37NT.cif", hotspots=HOTSPOTS)
    current = contacts(path, reference)
    conserved = sorted(baseline & current)
    moved = displacement(path, start)
    row.update(
        native_contact=closest_pair(path),
        kept_Ig1_contacts=conserved,
        baseline_Ig1_contacts=len(baseline),
        binder_pose_ca_rmsd_A=moved,
        pose_pass=bool(
            row["preliminary_geometry_pass"]
            and len(conserved) >= 0.7 * len(baseline)
            and moved <= max_rmsd
        ),
    )
    return row


def relax(pose, mode, pr):
    multisite = mode in MULTISITE_CASES
    cartesian = "backbone" in mode or multisite
    centers = (15, 32, 83) if multisite else (32,)
    ligand32 = pose.pdb_info().pdb2pose("B", 32)
    lys112 = pose.pdb_info().pdb2pose("A", 112)
    if (
        not ligand32
        or pose.residue(ligand32).name1() != "I"
        or not lys112
        or pose.residue(lys112).name1() != "K"
    ):
        raise ValueError("Expected binder Ile32 and wild-type receptor Lys112")
    mutations = (
        MULTISITE_CASES[mode] if multisite else {32: "GLY"} if mode.startswith("i32g") else {}
    )
    for number, identity in mutations.items():
        index = pose.pdb_info().pdb2pose("B", number)
        pr.rosetta.protocols.simple_moves.MutateResidue(index, identity).apply(pose)

    move = pr.rosetta.core.kinematics.MoveMap()
    move.set_bb(False)
    move.set_jump(False)
    move.set_chi(False)
    for idx in range(1, pose.total_residue() + 1):
        if pose.pdb_info().chain(idx) == "B":
            number = pose.pdb_info().number(idx)
            offset = min(abs(number - center) for center in centers)
            move.set_chi(idx, offset <= 6)
            move.set_bb(idx, cartesian and offset <= 3)
    move.set_chi(lys112, mode == "lys112_rotamer")
    rosetta = pr.rosetta
    score = pr.create_score_function("ref2015_cart" if cartesian else "ref2015")
    if cartesian:
        score.set_weight(rosetta.core.scoring.coordinate_constraint, 1.0)
        anchor = rosetta.core.id.AtomID(pose.residue(lys112).atom_index("CA"), lys112)
        for num in sorted({n for c in centers for n in range(c - 5, c + 6)}):
            idx = pose.pdb_info().pdb2pose("B", num)
            for atom in ("N", "CA", "C", "O"):
                aid = rosetta.core.id.AtomID(pose.residue(idx).atom_index(atom), idx)
                ideal = pose.xyz(aid)
                width = 0.6 if min(abs(num - c) for c in centers) <= 3 else 0.15
                func = rosetta.core.scoring.func.HarmonicFunc(0, width)
                pose.add_constraint(
                    rosetta.core.scoring.constraints.CoordinateConstraint(aid, anchor, ideal, func)
                )
    mover = rosetta.protocols.relax.FastRelax(score, 2)
    mover.cartesian(cartesian)
    mover.set_movemap(move)
    mover.set_movemap_disables_packing_of_fixed_chi_positions(True)
    mover.ramp_down_constraints(False)
    mover.max_iter(150)
    mover.apply(pose)


def screen(suite, start, baseline, cases):
    import pyrosetta as pr

    pr.init("-mute all -constant_seed -jran 112026")
    original = pr.pose_from_pdb(str(start))
    rows = {}
    for name in cases:
        print("Local test", name, flush=True)
        output = suite / f"{name}.pdb"
        pose = original.clone()
        relax(pose, name, pr)
        pose.dump_pdb(str(output))
        row = review(output, suite, start, baseline, max_rmsd=1.0)
        # A rotamer sensitivity test is scored against its own receptor as well as 37NT.
        print(name, row["native_contact"], row["pose_pass"], flush=True)
        rows[name] = row
        write_json(suite / "geometry.json", rows)


def check_refolds(suite, project, start, baseline):
    rows = json.loads((suite / "geometry.json").read_text())
    eligible = {name: row for name, row in rows.items() if row["pose_pass"]}
    # Refolding uses a receptor-only template, so repeat the original sequence only once.
    parents = {}
    for name, row in eligible.items():
        seq = row["binder_sequence"]
        if (
            seq not in parents
            or row["intact_protein_min_heavy_atom_distance_A"]
            > eligible[parents[seq]]["intact_protein_min_heavy_atom_distance_A"]
        ):
            parents[seq] = name
    results = {}
    for name in parents.values():
        path = suite / f"{name}.pdb"
        if file_hash(path) != eligible[name]["sha256"]:
            raise ValueError("Geometry candidate changed")
        output = suite / "refolds" / name
        output.mkdir(parents=True, exist_ok=True)
        reference = suite / "target.pdb"
        refold(
            project,
            reference,
            suite / "37NT.cif",
            read_design_chain(reference, "A")[0],
            read_design_chain(path, "B")[0],
            output,
            HOTSPOTS,
        )
        row = review(output / "refold.pdb", suite, start, baseline, max_rmsd=2.5)
        row["prediction"] = json.loads((output / "metrics.json").read_text())
        results[name] = row
        write_json(suite / "independent.json", results)
    write_json(suite / "independent.json", results)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["geometry", "multisite", "refold"])
    parser.add_argument("suite", type=Path)
    parser.add_argument("project", type=Path)
    args = parser.parse_args()
    start = args.suite / "mpnn8_parent_placed.pdb"
    baseline = {n for n in contacts(start, args.suite / "target.pdb") if 93 <= n <= 185}
    if args.stage != "refold":
        cases = GEOMETRY_CASES if args.stage == "geometry" else MULTISITE_CASES
        screen(args.suite, start, baseline, cases)
    else:
        check_refolds(args.suite, args.project, start, baseline)
