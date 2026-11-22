#!/usr/bin/env python3
"""Read an embedded lab-meeting snapshot and audit displayed Ig1 poses, without executing HTML."""

import argparse
import hashlib
import io
import json
import re
from pathlib import Path

import numpy as np
from Bio.PDB import PDBIO, MMCIFParser, PDBParser

from mertk.structural import rigid_fit, rmsd


def digest(data):
    return hashlib.sha256(data).hexdigest()


def embedded_json(html, name):
    matches = list(re.finditer(rf"const {name}\s*=\s*", html))
    if len(matches) != 1:
        raise ValueError(f"Expected one {name} assignment, found {len(matches)}")
    return json.JSONDecoder().raw_decode(html[matches[0].end() :])[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("html", type=Path)
    parser.add_argument(
        "--output", type=Path, default=Path("outputs/analysis/20261113/lab_reference")
    )
    args = parser.parse_args()
    raw = args.html.read_bytes()
    pages = embedded_json(raw.decode(), "PAGES")
    structure_page = next(p for p in pages if p["title"] == "MERTK structures")
    data = embedded_json(structure_page["html"], "DATA")
    reference_path = Path("data/references/37NT/37NT.cif")
    native = MMCIFParser(QUIET=True).get_structure("native", reference_path)[0]
    target = {r.id[1]: r for r in native["A"] if r.id[0] == " " and "CA" in r}
    glycan_atoms = [
        a
        for r in native.get_residues()
        if r.resname in {"NAG", "MAN", "BMA"}
        for a in r
        if a.element != "H"
    ]
    glycan_xyz = np.array([a.coord for a in glycan_atoms])
    attachment_atom = min(
        (a for a in glycan_atoms if a.parent.resname == "NAG" and a.name == "C1"),
        key=lambda a: np.linalg.norm(a.coord - target[114]["ND2"].coord),
    )
    species_path = Path("outputs/species_review_20261111/residue_mapping.json")
    species = {r["human_position"]: r for r in json.loads(species_path.read_text())}
    args.output.mkdir(parents=True, exist_ok=True)
    report = {
        "source_html": str(args.html.resolve()),
        "source_html_sha256": digest(raw),
        "structure_snapshot": data["snapshot"],
        "reference_sha256": {
            str(p): digest(p.read_bytes()) for p in [reference_path, species_path]
        },
        "reviewed_page": structure_page["title"],
        "n114_attachment_check": {
            "glycan_chain": attachment_atom.parent.parent.id,
            "glycan_residue": attachment_atom.parent.id[1],
            "ND2_to_NAG_C1_A": float(
                np.linalg.norm(attachment_atom.coord - target[114]["ND2"].coord)
            ),
        },
        "groups_reported": [
            {k: v for k, v in g.items() if k != "candidates"} for g in data["groups"]
        ],
        "controls_reported": [
            {
                "name": c["name"],
                "summary": c["summary"],
                "embedded_pose_count": len(c["candidates"]),
            }
            for c in data["controls"]
        ],
        "candidates": [],
        "limitations": [
            "Only the MERTK structure page is analyzed; this is a selected snapshot, not the complete campaign.",
            "Passing labels and metrics are reported by the snapshot; original filter/configuration files are unavailable.",
            "Displayed sequence variants share two trajectory hashes; no independent-attempt acceptance rate is inferred.",
            "Distances are post hoc rigid overlays onto one insect-cell glycosylated reference; no glycan sampling or energy minimization.",
            "Ig1 and Ig-pair fits use all matched CA atoms without outlier removal; high RMSD limits literal clash interpretation.",
            "Human/mouse sequence correspondence does not establish cross-species binding.",
            "Natural-ligand contact overlap uses one embedded predicted pose per control, not the reported three-seed ensemble.",
        ],
    }
    for group in data["groups"]:
        for c in group["candidates"]:
            structure = PDBParser(QUIET=True).get_structure(c["design"], io.StringIO(c["pdb"]))
            model = structure[0]
            predicted = {r.id[1]: r for r in model["A"] if r.id[0] == " " and "CA" in r}
            common = sorted(set(target) & set(predicted))
            mismatches = [n for n in common if target[n].resname != predicted[n].resname]
            if mismatches:
                raise ValueError(f"Target sequence mismatch at {mismatches}")
            binder = [a for r in model["B"] if r.id[0] == " " for a in r if a.element != "H"]
            binder_xyz = np.array([a.coord for a in binder])
            observed_contacts = []
            for residue in model["A"]:
                coordinates = np.array([a.coord for a in residue if a.element != "H"])
                if np.linalg.norm(coordinates[:, None] - binder_xyz[None, :], axis=-1).min() < 4.5:
                    observed_contacts.append(residue.id[1])
            if observed_contacts != c["interface"]:
                raise ValueError(
                    f"Displayed contact list disagrees with coordinates: {c['design']}"
                )
            record = {k: v for k, v in c.items() if k not in {"pdb", "contacts", "binderInterface"}}
            record["displayed_contacts_verified_at_4_5_A"] = True
            record["embedded_pdb_sha256"] = digest(c["pdb"].encode())
            record["size_group"] = group["size"]
            record["queued_hotspot_overlap"] = sorted(set(c["interface"]) & {122, 147, 148})
            record["human_mouse_contact_changes"] = [
                species[n]
                for n in c["interface"]
                if n in species and species[n]["human_aa"] != species[n]["mouse_aa"]
            ]
            record["ligand_contact_overlap"] = {
                control["name"]: sorted(
                    set(c["interface"]) & set(control["candidates"][0]["interface"])
                )
                for control in data["controls"]
                if control["name"] in {"GAS6", "Protein S"}
            }
            record["alignment_checks"] = {}
            for label, end in [("Ig1", 186), ("Ig_pair", 285)]:
                positions = [n for n in common if 93 <= n <= end]
                x = np.array([predicted[n]["CA"].coord for n in positions])
                y = np.array([target[n]["CA"].coord for n in positions])
                rotation, translation = rigid_fit(x, y)
                transformed = binder_xyz @ rotation + translation
                distances = np.linalg.norm(transformed[:, None] - glycan_xyz[None, :], axis=-1)
                bi, gi = np.unravel_index(distances.argmin(), distances.shape)
                ga = glycan_atoms[gi]
                record["alignment_checks"][label] = {
                    "ca_positions": positions,
                    "target_ca_rmsd_A": rmsd(x @ rotation + translation, y),
                    "min_binder_glycan_heavy_A": float(distances[bi, gi]),
                    "binder_atoms_within_2_5_A_of_glycan": int((distances.min(axis=1) < 2.5).sum()),
                    "closest_binder_atom": [binder[bi].parent.id[1], binder[bi].name],
                    "closest_glycan_atom": [
                        ga.parent.parent.id,
                        ga.parent.id[1],
                        ga.parent.resname,
                        ga.name,
                    ],
                    "rotation": rotation.tolist(),
                    "translation": translation.tolist(),
                }
                if label == "Ig1" and c["rank"] == "1":
                    aligned = structure.copy()
                    for atom in aligned.get_atoms():
                        atom.coord = atom.coord @ rotation + translation
                    output = args.output / f"lab_Ig1_l{c['length']}_rank1_aligned.pdb"
                    writer = PDBIO()
                    writer.set_structure(aligned)
                    writer.save(str(output))
                    record["visualization_copy"] = str(output)
            report["candidates"].append(record)
    report["displayed_candidate_count"] = len(report["candidates"])
    report["displayed_trajectory_hash_count"] = len({c["hash"] for c in report["candidates"]})
    (args.output / "reference_review.json").write_text(json.dumps(report, indent=2) + "\n")
    print(args.output / "reference_review.json")
    for c in report["candidates"]:
        print(
            c["length"],
            c["rank"],
            c["queued_hotspot_overlap"],
            {
                k: (
                    round(v["target_ca_rmsd_A"], 2),
                    round(v["min_binder_glycan_heavy_A"], 2),
                    v["closest_glycan_atom"],
                )
                for k, v in c["alignment_checks"].items()
            },
            "mouse changes",
            c["human_mouse_contact_changes"],
        )


if __name__ == "__main__":
    main()
