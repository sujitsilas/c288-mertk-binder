"""Exploratory MERTK stalk geometry; no structure score establishes shedding protection."""

import copy
import json
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFIO, PDBIO, MMCIFParser, PDBParser, Select
from Bio.SeqUtils import seq1

from mertk.glycan_guard import MODE, SUGARS, frames
from mertk.prediction_review import nearest_points
from mertk.protector_screen import bond_distances
from mertk.provenance import file_hash
from mertk.structural import rigid_fit, rmsd

TARGET_RANGE = (386, 526)
STALK_RANGE = (485, 498)
HOTSPOTS = [488, 490, 491]


class Crop(Select):
    def __init__(self, low, high):
        self.low, self.high = low, high

    def accept_chain(self, chain):
        return chain.id == "A"

    def accept_residue(self, residue):
        return residue.id[0] == " " and self.low <= residue.id[1] <= self.high


def heavy(chain):
    return [a for r in chain if r.id[0] == " " for a in r if a.element not in {"H", "D"}]


def sequence(chain):
    return "".join(seq1(r.resname) for r in chain if r.id[0] == " ")


def rotate_psi(model, residue, degrees):
    """Perturb one backbone torsion; this is sensitivity sampling, not an equilibrium ensemble."""
    model = copy.deepcopy(model)
    chain = model["A"]
    origin = chain[residue]["C"].coord.copy()
    axis = origin - chain[residue]["CA"].coord
    axis /= np.linalg.norm(axis)
    theta = np.deg2rad(degrees)
    for r in chain:
        for atom in r:
            if r.id[1] > residue or (r.id[1] == residue and atom.name in {"O", "OXT"}):
                x = atom.coord - origin
                atom.coord = (
                    x * np.cos(theta)
                    + np.cross(axis, x) * np.sin(theta)
                    + axis * np.dot(axis, x) * (1 - np.cos(theta))
                    + origin
                )
    return model


def glycan_constraints(target_path, native_cif):
    target = PDBParser(QUIET=True).get_structure("model", target_path)[0]["A"]
    numbers = [r.id[1] for r in target if r.id[0] == " "]
    if numbers != list(range(TARGET_RANGE[0], TARGET_RANGE[1] + 1)):
        raise ValueError("Protector requires the declared contiguous A386–A526 context")
    native = MMCIFParser(QUIET=True).get_structure("native", native_cif)[0]
    local, atom_sites, sites, attachments = [], [], [], []
    for chain in native:
        sugars = [r for r in chain if r.resname in SUGARS]
        if not sugars:
            continue
        candidates = sorted(
            (float(np.linalg.norm(s["C1"].coord - r["ND2"].coord)), r.id[1])
            for s in sugars
            if s.resname == "NAG" and "C1" in s
            for r in native["A"]
            if r.resname == "ASN" and "ND2" in r
        )
        if not candidates or candidates[0][0] > 2.2:
            raise ValueError("Unresolved carbohydrate attachment")
        distance, n = candidates[0]
        if n not in numbers:
            continue
        triad = np.array([native["A"][n][a].coord for a in ["N", "CA", "C"]])
        xyz = np.array([a.coord for r in sugars for a in r if a.element not in {"H", "D"}])
        local.extend(((xyz - triad[1]) @ frames(triad).T).tolist())
        atom_sites.extend([len(sites)] * len(xyz))
        sites.append(n)
        attachments.append(
            dict(
                native_asn=n,
                glycan_chain=chain.id,
                covalent_distance_A=distance,
                heavy_atoms=len(xyz),
            )
        )
    residues = list(target)
    sequons = [
        r.id[1]
        for i, r in enumerate(residues[:-2])
        if r.resname == "ASN"
        and residues[i + 1].resname != "PRO"
        and residues[i + 2].resname in {"SER", "THR"}
    ]
    if sorted(sites) != [389, 395] or sequons != [389, 395, 442, 454]:
        raise ValueError("Unexpected protector glycan or sequon mapping")
    return dict(
        mode=MODE,
        target_sha256=file_hash(target_path),
        native_cif_sha256=file_hash(native_cif),
        target_length=len(numbers),
        native_numbers=numbers,
        attachments=attachments,
        anchor_indices=[numbers.index(n) for n in sites],
        sequon_indices=[numbers.index(n) for n in sequons],
        glycan_local_xyz=local,
        glycan_atom_sites=atom_sites,
        heavy_atom_clearance_A=3.5,
        sequon_ca_exclusion_A=15.0,
        weights={"glycan_volume": 10.0, "sequon_volume": 2.0, "membrane_volume": 3.0},
        membrane={
            "origin_index": numbers.index(506),
            "outer_indices": [numbers.index(n) for n in range(506, 513)],
            "inner_indices": [numbers.index(n) for n in range(520, 527)],
            "extracellular_margin_A": 5.0,
        },
        limitations="Modeled stalk and membrane orientation; locally transported observed glycans, "
        "unresolved sequons excluded conservatively. No membrane or glycan cofolding.",
    )


def transported_glycans(chain, spec):
    residues = [chain[n] for n in spec["native_numbers"]]
    triads = np.array(
        [[residues[i][a].coord for a in ["N", "CA", "C"]] for i in spec["anchor_indices"]]
    )
    ids = np.asarray(spec["glycan_atom_sites"])
    return (
        np.einsum("ai,aij->aj", np.asarray(spec["glycan_local_xyz"]), frames(triads)[ids])
        + triads[ids, 1]
    )


def membrane_clearance(xyz, chain):
    outer = np.mean([chain[n]["CA"].coord for n in range(506, 513)], axis=0)
    inner = np.mean([chain[n]["CA"].coord for n in range(520, 527)], axis=0)
    axis = (inner - outer) / np.linalg.norm(inner - outer)
    return float(-np.max((xyz - chain[506]["CA"].coord) @ axis) - 5.0)


def prepare_targets(source, native_cif, canonical, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    original = PDBParser(QUIET=True).get_structure("af", source)[0]
    if sequence(original["A"]) != canonical:
        raise ValueError("AlphaFold model sequence differs from canonical MERTK")
    reports = []
    for name, angle in [("base", 0), ("psi490_plus25", 25), ("psi490_minus25", -25)]:
        model = rotate_psi(original, 490, angle)
        writer = PDBIO()
        writer.set_structure(model)
        target = out / f"{name}_target.pdb"
        writer.save(str(target), Crop(*TARGET_RANGE))
        writer.save(str(out / f"{name}_full_context.pdb"), Crop(93, 526))
        constraint = glycan_constraints(target, native_cif)
        (out / f"{name}_constraints.json").write_text(json.dumps(constraint, indent=2) + "\n")
        chain = model["A"]
        ca = np.array([chain[n]["CA"].coord for n in range(386, 527)])
        steps = np.linalg.norm(np.diff(ca, axis=0), axis=1)
        if (
            not np.all((steps > 2.7) & (steps < 4.5))
            or not 1.2 < chain[490]["C"] - chain[491]["N"] < 1.5
        ):
            raise ValueError("Conformer has broken peptide geometry")
        bond = np.mean([chain[490]["C"].coord, chain[491]["N"].coord], axis=0)
        gly = transported_glycans(chain, constraint)
        reports.append(
            {
                "conformer": name,
                "psi490_delta_deg": angle,
                "cleavage_bond_length_A": float(chain[490]["C"] - chain[491]["N"]),
                "stalk_plddt": [chain[n]["CA"].bfactor for n in range(485, 499)],
                "bond_to_observed_glycan_min_A": float(np.linalg.norm(gly - bond, axis=1).min()),
                "bond_membrane_clearance_with_5A_margin": membrane_clearance(bond[None], chain),
                "target_sha256": file_hash(target),
            }
        )
    # Model-only visualization reference, never labeled experimental.
    visual = copy.deepcopy(original)
    for chain in list(visual):
        if chain.id != "A":
            visual.detach_child(chain.id)
    for r in list(visual["A"]):
        if not 93 <= r.id[1] <= 526:
            visual["A"].detach_child(r.id)
    writer = MMCIFIO()
    writer.set_structure(visual)
    writer.save(str(out / "predicted_context.cif"))
    report = {
        "target": "Human MERTK Q12866",
        "hypothesis": "Monovalent occupancy of the P490–S491 cleavage region",
        "hotspots": HOTSPOTS,
        "target_range": TARGET_RANGE,
        "binding_region": STALK_RANGE,
        "experimental_structure_ends_at": 484,
        "source_sha256": file_hash(source),
        "conformers": reports,
        "generation_conformers": ["base", "psi490_plus25"],
        "sensitivity_only_conformer": "psi490_minus25",
        "limitations": "Low-confidence stalk. Psi perturbations test geometric sensitivity and are not "
        "a physical conformational ensemble. Binding or AF3 confidence does not prove cleavage inhibition.",
    }
    (out / "triage.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def evaluate_pose(path, campaign, suite):
    model = PDBParser(QUIET=True).get_structure("prediction", path)[0]
    reference = PDBParser(QUIET=True).get_structure("reference", campaign / "inputs/target.pdb")[0][
        "A"
    ]
    target = next((c for c in model if sequence(c) == sequence(reference)), None)
    if target is None or len(list(model)) != 2:
        raise ValueError("Expected exact target sequence and one binder")
    binder = next(c for c in model if c.id != target.id)
    mapped = {n: r for n, r in zip(range(386, 527), target)}
    xyz = np.array([a.coord for a in heavy(binder)])
    contacts = [
        n
        for n, r in mapped.items()
        if nearest_points(xyz, np.array([a.coord for a in r if a.element != "H"]))[0].min() <= 4.5
    ]
    stalk = [n for n in contacts if 485 <= n <= 498]
    distances = bond_distances(xyz, mapped)
    reasons = []
    if max(distances) > 5 or len(stalk) < 3 or len(stalk) < 0.5 * max(1, len(contacts)):
        reasons.append("Insufficient cleavage-site contact or predominantly off-stalk interface")
    if any(n >= 506 for n in contacts):
        reasons.append("Transmembrane contacts")
    checks = []
    spec = json.loads((campaign / "inputs/glycan_constraints.json").read_text())
    native = MMCIFParser(QUIET=True).get_structure("native", suite / "references/37NT.cif")[0]["A"]
    for name in ["base", "psi490_plus25", "psi490_minus25"]:
        context = PDBParser(QUIET=True).get_structure(
            "context", suite / f"references/{name}_full_context.pdb"
        )[0]["A"]
        ns = range(485, 499)
        x = np.array([mapped[n]["CA"].coord for n in ns])
        y = np.array([context[n]["CA"].coord for n in ns])
        rotation, translation = rigid_fit(x, y)
        moved = xyz @ rotation + translation
        pd, _ = nearest_points(moved, np.array([a.coord for a in heavy(context)]))
        gd, _ = nearest_points(moved, transported_glycans(context, spec))
        nd = range(386, 480)
        nr, nt = rigid_fit(
            np.array([native[n]["CA"].coord for n in nd]),
            np.array([context[n]["CA"].coord for n in nd]),
        )
        experimental = np.array([a.coord for a in heavy(native)]) @ nr + nt
        ed, _ = nearest_points(moved, experimental)
        sd, _ = nearest_points(
            moved, np.array([context[n]["CA"].coord for n in [389, 395, 442, 454]])
        )
        row = dict(
            conformer=name,
            stalk_ca_rmsd_A=rmsd(x @ rotation + translation, y),
            model_protein_min_A=float(pd.min()),
            restored_protein_min_A=float(ed.min()),
            glycan_min_A=float(gd.min()),
            sequon_ca_min_A=float(sd.min()),
            membrane_clearance_A=membrane_clearance(moved, context),
        )
        row["pass"] = bool(
            row["stalk_ca_rmsd_A"] <= 2.5
            and pd.min() >= 2
            and ed.min() >= 2
            and gd.min() >= 2.5
            and sd.min() >= 10
            and row["membrane_clearance_A"] >= 0
        )
        checks.append(row)
    if not all(c["pass"] for c in checks):
        reasons.append("Failed one or more conformer/context clearance checks")
    return dict(
        path=str(path),
        sha256=file_hash(path),
        binder_sequence=sequence(binder),
        contacts=contacts,
        stalk_contacts=stalk,
        cleavage_atom_min_distances_A=distances,
        conformers=checks,
        rejection_reasons=reasons,
        preliminary_model_pass=not reasons,
        qualified=False,
        cleavage_protection_predicted=False,
    )
