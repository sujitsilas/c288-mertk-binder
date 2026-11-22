#!/usr/bin/env python3
"""Test native and alternative MerTK rotamers around a Boltz-recovered Ig1 pose."""

import argparse
import copy
import json
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser

from mertk.prediction_review import SUGARS, native_numbering, pdb_chain, sequence
from mertk.provenance import file_hash
from mertk.structural import rigid_fit, rmsd

LYSINES = (109, 112)
ANCHORS = (122, 148)
BACKBONE = ("N", "CA", "C", "O")
SEEDS = (202, 303)
VARIANTS = {
    "e80d": {80: "ASP"},
    "e80q": {80: "GLN"},
    "f83v": {83: "VAL"},
    "e80d_f83v": {80: "ASP", 83: "VAL"},
}


def coords(atoms):
    return np.array([atom.coord for atom in atoms], dtype=float)


def heavy(atoms):
    return [atom for atom in atoms if atom.element not in {"H", "D"}]


def proteins(chain):
    return [r for r in chain if r.id[0] == " "]


def distance(left, right):
    squared = (left**2).sum(1)[:, None] + (right**2).sum(1) - 2 * left @ right.T
    return np.maximum(0, squared) ** 0.5


def save_json(path, value):
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    tmp.replace(path)


def prepare(source, out):
    """Place the Boltz binder and local lysine alternatives in the 37NT frame."""
    out.mkdir(parents=True, exist_ok=False)
    manifest = json.loads((source / "manifest.json").read_text())
    name = "01_ig1_i32g_w83f"
    case = next(c for c in manifest["cases"] if c["id"] == name)
    rows = json.loads((source / "assessment.json").read_text())["cases"]
    row = next(r for r in rows if r["id"] == name)
    raw_path = source / row["structure"]
    canonical_path = source / "references/canonical.json"
    references = source / "references/37NT.cif"
    if file_hash(raw_path) != row["source_sha256"]:
        raise ValueError("Scored Boltz structure has changed")
    if file_hash(canonical_path) != manifest["canonical_source_sha256"]:
        raise ValueError("Canonical sequence has changed")
    canonical = json.loads(canonical_path.read_text())["sequence"]["value"]
    parser = MMCIFParser(QUIET=True)
    native = parser.get_structure("37NT", references)[0]["A"]
    raw = parser.get_structure("Boltz", raw_path)[0]
    if (sequence(raw["A"]), sequence(raw["B"])) != (
        manifest["receptor_sequence"],
        case["sequence"],
    ):
        raise ValueError("Boltz model does not match the frozen chain sequences")
    mapping = native_numbering(raw["A"], canonical)
    predicted = {mapping[r.id]: r for r in proteins(raw["A"])}
    original = {r.id[1]: r for r in proteins(native)}
    positions = range(93, 186)
    transform = rigid_fit(
        coords(predicted[n]["CA"] for n in positions), coords(original[n]["CA"] for n in positions)
    )
    partner = pdb_chain(raw["B"], transform=transform)
    modeled = copy.deepcopy(native)
    for n in LYSINES:
        a, b = predicted[n], modeled[n]
        anchors = ("N", "CA", "C")
        rotation, translation = rigid_fit(
            coords(a[k] for k in anchors), coords(b[k] for k in anchors)
        )
        for atom in b:
            if atom.name not in BACKBONE and atom.element != "H":
                atom.coord = a[atom.name].coord @ rotation + translation
    for label, receptor in [("crystal", native), ("boltz_lysines", modeled)]:
        receptor_block = pdb_chain(receptor).rsplit("END", 1)[0]
        (out / f"{label}.pdb").write_text(receptor_block + partner)
    (out / "37NT.cif").write_bytes(references.read_bytes())
    (out / "canonical.json").write_bytes(canonical_path.read_bytes())
    save_json(
        out / "inputs.json",
        {
            "case": name,
            "source_sha256": row["source_sha256"],
            "expected_binder": case["sequence"],
            "expected_receptor": sequence(native),
            "starting_contacts": case["reference_contact_residues_4p5A"],
            "templates": ["crystal", "boltz_lysines"],
            "variant_sites": VARIANTS,
            "seeds": list(SEEDS),
            "inputs": {f.name: file_hash(f) for f in sorted(out.iterdir())},
        },
    )


def geometry(model, native, start, starting_contacts, binder_limit=0.05):
    """Compare the saved complex, fixed-frame crystal context and observed glycans."""
    receptor = heavy(a for r in proteins(model["A"]) for a in r)
    ligand = heavy(model["B"].get_atoms())
    ref = heavy(a for r in proteins(native["A"]) for a in r)
    glycans = heavy(a for r in native.get_residues() if r.resname in SUGARS for a in r)
    receptor_ids = np.array([a.parent.id[1] for a in receptor])
    xyz = coords(ligand)
    distances = distance(coords(receptor), xyz)
    d_native = distance(coords(ref), xyz)
    glycan = float(distance(coords(glycans), xyz).min())
    contacts = receptor_ids[distances.min(axis=1) <= 4.5].tolist()
    ig1 = sorted(set(contacts).intersection(range(93, 186)))
    retained = sorted(set(starting_contacts).intersection(ig1))
    hotspots = float(distances[np.isin(receptor_ids, ANCHORS)].min())
    receptor_ca = rmsd(
        coords(r["CA"] for r in proteins(model["A"])),
        coords(r["CA"] for r in proteins(native["A"])),
    )
    binder_ca = rmsd(
        coords(r["CA"] for r in proteins(model["B"])), coords(r["CA"] for r in proteins(start["B"]))
    )
    if receptor_ca > 0.05:
        raise ValueError("Repacking altered the fixed receptor backbone frame")
    return {
        "closest_protein_A": float(distances.min()),
        "pairs_below_2A": int((distances < 2).sum()),
        "fixed_37NT_min_A": float(d_native.min()),
        "fixed_37NT_pairs_below_2A": int((d_native < 2).sum()),
        "observed_glycan_min_A": glycan,
        "retained_Ig1_contacts": retained,
        "Ig1_contacts": ig1,
        "starting_contact_count": len(starting_contacts),
        "anchor_min_A": hotspots,
        "receptor_ca_rmsd_A": receptor_ca,
        "binder_ca_rmsd_A": binder_ca,
        "geometry_pass": bool(
            distances.min() >= 2
            and glycan >= 2.5
            and len(retained) >= 0.7 * len(starting_contacts)
            and hotspots <= 5
            and binder_ca <= binder_limit
        ),
    }


def relax(pose, sites, score, pr, redock=False):
    rosetta = pr.rosetta
    move = rosetta.core.kinematics.MoveMap()
    move.set_bb(False)
    move.set_jump(False)
    move.set_chi(False)
    for index in sites:
        move.set_chi(index, True)
    if redock:
        score = score.clone()
        score.set_weight(rosetta.core.scoring.coordinate_constraint, 0.6)
        move.set_jump(binder_jump_index(pose), True)
        info = pose.pdb_info()
        anchor = info.pdb2pose("A", 109)
        anchor_atom = rosetta.core.id.AtomID(pose.residue(anchor).atom_index("CA"), anchor)
        for n in range(1, 126, 5):
            i = info.pdb2pose("B", n)
            atom = rosetta.core.id.AtomID(pose.residue(i).atom_index("CA"), i)
            restraint = rosetta.core.scoring.constraints.CoordinateConstraint(
                atom, anchor_atom, pose.xyz(atom), rosetta.core.scoring.func.HarmonicFunc(0, 0.6)
            )
            pose.add_constraint(restraint)
    mover = rosetta.protocols.relax.FastRelax(score, 1)
    mover.set_movemap(move)
    mover.set_movemap_disables_packing_of_fixed_chi_positions(True)
    mover.ramp_down_constraints(False)
    mover.max_iter(120)
    mover.apply(pose)


def binder_jump_index(pose):
    tree = pose.fold_tree()
    info = pose.pdb_info()
    (binder_jump,) = [
        i for i in range(1, tree.num_jump() + 1) if info.chain(tree.jump_edge(i).stop()) == "B"
    ]
    return binder_jump


def interface_terms(pose, score, pr):
    terms = pr.rosetta.core.scoring
    total = score(pose)
    repulsion = pose.energies().total_energies()[terms.fa_rep]
    separated = pose.clone()
    info = separated.pdb_info()
    binder_jump = binder_jump_index(separated)
    jump = separated.jump(binder_jump)
    jump.set_translation(jump.get_translation() + pr.rosetta.numeric.xyzVector_double_t(1000, 0, 0))
    separated.set_jump(binder_jump, jump)
    binder_marker = info.pdb2pose("B", 1)
    native_marker = info.pdb2pose("A", 109)
    sugar_marker = info.pdb2pose("A", 601)
    if (
        separated.residue(binder_marker).xyz("CA").distance(pose.residue(binder_marker).xyz("CA"))
        < 500
        or separated.residue(native_marker)
        .xyz("CA")
        .distance(pose.residue(native_marker).xyz("CA"))
        > 0.01
        or (
            sugar_marker
            and separated.residue(sugar_marker).xyz(1).distance(pose.residue(sugar_marker).xyz(1))
            > 0.01
        )
    ):
        raise ValueError("Rigid separation moved the wrong chain")
    separated_total = score(separated)
    separated_repulsion = separated.energies().total_energies()[terms.fa_rep]
    return {
        "interface_delta_REU": total - separated_total,
        "interface_fa_rep_unweighted": repulsion - separated_repulsion,
        "lys_fa_dun_unweighted": {
            n: pose.energies().residue_total_energies(pose.pdb_info().pdb2pose("A", n))[
                terms.fa_dun
            ]
            for n in LYSINES
        },
    }


def rescore(out):
    """Measure the interface against a rigidly separated copy of the same complex."""
    import pyrosetta as pr

    pr.init("-mute all -no_optH")
    score = pr.get_fa_scorefxn()
    rows = json.loads((out / "repack.json").read_text())["cases"]
    results = []
    for row in rows:
        path = out / row["file"]
        if file_hash(path) != row["sha256"]:
            raise ValueError(f"Frozen Rosetta output changed: {path}")
        pose = pr.pose_from_pdb(str(path))
        energy = interface_terms(pose, score, pr)
        results.append({**row, **energy})
        print(row["id"], energy, flush=True)
    save_json(
        out / "repack_energy.json",
        {
            "method": "REF2015 bound minus rigidly separated complex, without repacking",
            "cases": results,
        },
    )


def polish(out, redock=False):
    """Resolve the residual His26 and Phe37 contacts with the pose frame fixed."""
    import pyrosetta as pr

    pr.init("-mute all -ex1 -ex2 -use_input_sc -constant_seed -jran 202")
    metadata = json.loads((out / "inputs.json").read_text())
    rows = json.loads((out / "repack.json").read_text())["cases"]
    native = MMCIFParser(QUIET=True).get_structure("37NT", out / "37NT.cif")[0]
    pdb = PDBParser(QUIET=True)
    score = pr.get_fa_scorefxn()
    results = []
    candidates = {"binder_only", "both_sides", "e80q", "f83v"}
    for parent in rows:
        if (
            parent["start"] != "crystal"
            or parent["seed"] != 202
            or parent["protocol"] not in candidates
        ):
            continue
        file = out / parent["file"]
        if file_hash(file) != parent["sha256"]:
            raise ValueError(f"Frozen repacked pose changed: {file}")
        if redock:
            previous = json.loads((out / "polish.json").read_text())["cases"]
            polished = next(row for row in previous if row["id"] == parent["id"])
            file = out / polished["file"]
            if file_hash(file) != polished["sha256"]:
                raise ValueError(f"Frozen polished pose changed: {file}")
        source = pdb.get_structure("start", file)[0]
        pose = pr.pose_from_pdb(str(file))
        info = pose.pdb_info()
        residues = [
            *[info.pdb2pose("B", n) for n in sorted({*parent["binder_repack_sites"], 26, 37})],
            *[info.pdb2pose("A", n) for n in parent["receptor_repack_sites"]],
        ]
        initial = interface_terms(pose, score, pr)
        relax(pose, residues, score, pr, redock=redock)
        stage = "redocked" if redock else "polished"
        polished_path = out / "repacked" / f"{parent['id']}_{stage}.pdb"
        pose.dump_pdb(str(polished_path))
        polished = pdb.get_structure("polished", polished_path)[0]
        hh_min = {}
        for label, model in [("before", source), ("after", polished)]:
            a, b = [
                coords(atom for atom in model[c].get_atoms() if atom.element == "H") for c in "AB"
            ]
            hh_min[label] = float(distance(a, b).min())
        row = {
            "id": parent["id"],
            "file": str(polished_path.relative_to(out)),
            "sha256": file_hash(polished_path),
            "hydrogen_pair_min_A": hh_min,
            "interface_before": initial,
            "interface_after": interface_terms(pose, score, pr),
            **geometry(
                polished, native, source, metadata["starting_contacts"], 1.0 if redock else 0.05
            ),
        }
        print(parent["id"], hh_min, row["closest_protein_A"], row["interface_after"], flush=True)
        results.append(row)
        save_json(out / ("redock.json" if redock else "polish.json"), {"cases": results})


def run(out):
    import pyrosetta as pr

    meta = json.loads((out / "inputs.json").read_text())
    for name, expected_hash in meta["inputs"].items():
        if file_hash(out / name) != expected_hash:
            raise ValueError(f"Frozen input changed: {name}")
    pr.init("-mute all -ex1 -ex2 -use_input_sc -constant_seed -jran 202")
    native = MMCIFParser(QUIET=True).get_structure("37NT", out / "37NT.cif")[0]
    parser = PDBParser(QUIET=True)
    score = pr.get_fa_scorefxn()
    results = []
    results_dir = out / "repacked"
    results_dir.mkdir(exist_ok=True)
    cases = [
        ("crystal", "binder_only", {}),
        ("crystal", "both_sides", {}),
        ("boltz_lysines", "both_sides", {}),
        *[("crystal", label, mutations) for label, mutations in VARIANTS.items()],
    ]
    for seed in SEEDS:
        for start_name, label, mutations in cases:
            start_path = out / f"{start_name}.pdb"
            start = parser.get_structure(start_name, start_path)[0]
            pose = pr.pose_from_pdb(str(start_path))
            info = pose.pdb_info()
            if (
                sequence(start["B"]) != meta["expected_binder"]
                or sequence(start["A"]) != meta["expected_receptor"]
            ):
                raise ValueError("Prepared complexes have an unexpected sequence")
            centers = [start["B"][i] for i in [80, 83]]
            region = coords(heavy(a for residue in centers for a in residue))
            binder_sites = [
                info.pdb2pose("B", res.id[1])
                for res in proteins(start["B"])
                if distance(coords(heavy(res)), region).min() <= 5.0
            ]
            receptor_sites = (
                [info.pdb2pose("A", n) for n in LYSINES] if label == "both_sides" else []
            )
            pr.rosetta.numeric.random.rg().set_seed(seed)
            for n, aa in mutations.items():
                index = info.pdb2pose("B", n)
                pr.rosetta.protocols.simple_moves.MutateResidue(index, aa).apply(pose)
            relax(pose, receptor_sites + binder_sites, score, pr)
            name = f"{start_name}_{label}_{seed}"
            path = results_dir / f"{name}.pdb"
            pose.dump_pdb(str(path))
            repacked = parser.get_structure(name, path)[0]
            if sequence(repacked["A"]) != meta["expected_receptor"]:
                raise ValueError("Repacked receptor sequence changed")
            row = {
                "id": name,
                "start": start_name,
                "protocol": label,
                "seed": seed,
                "file": str(path.relative_to(out)),
                "sha256": file_hash(path),
                "binder_sequence": sequence(repacked["B"]),
                "binder_repack_sites": [info.number(i) for i in binder_sites],
                "receptor_repack_sites": [info.number(i) for i in receptor_sites],
                **geometry(repacked, native, start, meta["starting_contacts"]),
                **interface_terms(pose, score, pr),
            }
            results.append(row)
            save_json(out / "repack.json", {"cases": results})
            print(
                name,
                row["geometry_pass"],
                round(row["closest_protein_A"], 2),
                round(row["fixed_37NT_min_A"], 2),
                len(row["retained_Ig1_contacts"]),
                round(row["interface_delta_REU"], 2),
                flush=True,
            )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["prepare", "run", "rescore", "polish", "redock"])
    parser.add_argument("suite", type=Path)
    parser.add_argument("--source", type=Path)
    args = parser.parse_args()
    if args.stage == "prepare":
        prepare(args.source, args.suite)
    elif args.stage == "run":
        run(args.suite)
    elif args.stage == "rescore":
        rescore(args.suite)
    else:
        polish(args.suite, redock=args.stage == "redock")


if __name__ == "__main__":
    main()
