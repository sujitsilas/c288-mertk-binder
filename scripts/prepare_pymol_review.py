#!/usr/bin/env python3
"""Assemble an auditable local structure review; proposed patches are not validated hotspots."""

import argparse
import csv
import json
import shutil
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser

from mertk.provenance import file_hash
from mertk.structural import rigid_fit, rmsd
from mertk.structures import read_design_chain


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    raw = out / "structures"
    raw.mkdir(exist_ok=True)
    source = Path("data/references/pymol_review")
    species = Path("outputs/species_review_20261111")
    campaign = Path("outputs/hpc/pilot_20261110_01")
    downloads = json.loads((source / "downloads.json").read_text())
    entries = [
        {
            "file": "37NT.cif",
            "url": "https://files.rcsb.org/download/37NT.cif",
            "type": "experimental human MERTK",
            "sha256": file_hash("data/references/37NT/37NT.cif"),
        }
    ]
    shutil.copy2("data/references/37NT/37NT.cif", raw / "37NT.cif")
    for entry in downloads:
        assert file_hash(source / entry["file"]) == entry["sha256"]
        shutil.copy2(source / entry["file"], raw / entry["file"])
        entries.append(entry)
    h = json.loads((species / "Q12866.json").read_text())["sequence"]["value"]
    mouse = {
        r["human_position"]: r
        for r in json.loads((species / "residue_mapping.json").read_text())
        if r["human_position"]
    }
    tam = json.loads((source / "tam_mapping.json").read_text())
    sasa = json.loads((source / "mertk_sasa.json").read_text())
    model = MMCIFParser(QUIET=True).get_structure("mertk", raw / "37NT.cif")[0]
    chain = model["A"]
    glycans = np.array(
        [
            a.coord
            for a in model.get_atoms()
            if a.parent.resname in {"NAG", "MAN", "BMA"} and a.element != "H"
        ]
    )
    sequons = [
        i + 1 for i in range(92, 482) if h[i] == "N" and h[i + 1] != "P" and h[i + 2] in "ST"
    ]
    patches = [
        {
            "id": "A_FNIIIa",
            "anchors": [365, 367, 369],
            "center": 367,
            "rationale": "Exposed FNIIIa flank, outside the Ig ligand-binding module; anchors conserved in mouse and differ at sequence-aligned positions in human AXL/TYRO3.",
            "limitations": "Nearby N316/N354 glycans may shield a larger binder footprint. No affinity or agonism evidence.",
        },
        {
            "id": "B_FNIIIb",
            "anchors": [460, 471, 472],
            "center": 472,
            "rationale": "Exposed alternative FNIIIb flank; anchors conserved in mouse. Weaker sequence-based TAM discrimination than patch A.",
            "limitations": "FNIIIb can participate in receptor contacts; folded-domain location does not demonstrate steric protection of the unresolved stalk.",
        },
    ]
    rows = []
    for patch in patches:
        patch["anchor_labels"] = [f"{h[n - 1]}{n}" for n in patch["anchors"]]
        for n in patch["anchors"]:
            atoms = np.array([a.coord for a in chain[n] if a.element != "H"])
            rows.append(
                {
                    "patch": patch["id"],
                    "human_residue": f"{h[n - 1]}{n}",
                    "human_position": n,
                    "mouse_position": mouse[n]["mouse_position"],
                    "mouse_amino_acid": mouse[n]["mouse_aa"],
                    "sasa_A2": round(sasa[str(n)]["sasa"], 2),
                    "nearest_observed_glycan_heavy_atom_A": round(
                        float(np.linalg.norm(atoms[:, None] - glycans[None, :], axis=-1).min()), 2
                    ),
                    "AXL_aligned_position": tam["P30530"][str(n)][0],
                    "AXL_amino_acid": tam["P30530"][str(n)][1],
                    "TYRO3_aligned_position": tam["Q06418"][str(n)][0],
                    "TYRO3_amino_acid": tam["Q06418"][str(n)][1],
                }
            )
    with (out / "candidate_patches.csv").open("w") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    reference, mapping = read_design_chain(campaign / "inputs/target.pdb", "A")
    numbers = np.array([r["author_number"] for r in mapping])
    diagnostics = []
    for f in sorted(campaign.glob("rfd4_proteina/inference_samples/*/*_generated.cif")):
        shutil.copy2(f, raw / f.name)
        entries.append(
            {
                "file": f.name,
                "type": "unvalidated RFD4 pilot output",
                "sha256": file_hash(f),
                "source": str(f),
            }
        )
        target, _ = read_design_chain(f, "A")
        assert target.sequence == reference.sequence
        rot, trans = rigid_fit(target.backbone[:, 1], reference.backbone[:, 1])
        row = {
            "file": f.name,
            "target_ca_rmsd": rmsd(target.backbone[:, 1] @ rot + trans, reference.backbone[:, 1]),
            "domains": {},
        }
        for domain, lo, hi in [
            ("Ig1", 93, 186),
            ("Ig2", 197, 273),
            ("FNIIIa", 286, 381),
            ("FNIIIb", 386, 484),
        ]:
            mask = (numbers >= lo) & (numbers <= hi)
            a, b = target.backbone[mask, 1], reference.backbone[mask, 1]
            r, t = rigid_fit(a, b)
            row["domains"][domain] = rmsd(a @ r + t, b)
        distances = np.linalg.norm(np.diff(target.backbone[:, 1], axis=0), axis=1)
        row["adjacent_ca_distance_range_A"] = [float(distances.min()), float(distances.max())]
        row["adjacent_ca_outliers_outside_3_to_4p5_A"] = [
            {"human_positions": [int(numbers[i]), int(numbers[i + 1])], "distance_A": float(d)}
            for i, d in enumerate(distances)
            if d < 3 or d > 4.5
        ]
        row["rotation"] = rot.tolist()
        row["translation"] = trans.tolist()
        diagnostics.append(row)
    reference_steps = np.linalg.norm(np.diff(reference.backbone[:, 1], axis=0), axis=1)
    report = {
        "date": "2026-11-11",
        "status": "Exploratory surface patches only; not submitted as design hotspots",
        "sources": entries,
        "patches": patches,
        "measurements": rows,
        "potential_N_glycosylation_sequons": sequons,
        "human_mouse_differences": [n for n in range(93, 485) if mouse[n]["mouse_aa"] != h[n - 1]],
        "experimental_glycan_sites": [114, 170, 207, 215, 234, 294, 316, 354, 389, 395],
        "crystal_contact_caution_residues": [444, 449],
        "methods": [
            "PyMOL solvent-accessible surface area; 1.4 A probe, dot_density=3, observed protein/glycan structure without water.",
            "Human-mouse canonical ectodomain and domain-specific human TAM sequence alignments: BLOSUM62, gap open -10, extension -0.5. Paralogue correspondence is approximate and not a selectivity prediction.",
            "Observed glycans are incomplete insect-cell glycans. All N-X-S/T sequons are marked, not assumed occupied.",
            "RFD4 domains aligned independently by sequence-matched CA Kabsch fits. Adjacent CA flags are geometric diagnostics, not an experimental quality threshold.",
        ],
        "rfd4_diagnostics": diagnostics,
        "experimental_reference_adjacent_ca_distance_range_A": [
            float(reference_steps.min()),
            float(reference_steps.max()),
        ],
        "experimental_reference_pairs_outside_3_to_4p5_A": int(
            ((reference_steps < 3) | (reference_steps > 4.5)).sum()
        ),
    }
    (out / "review.json").write_text(json.dumps(report, indent=2) + "\n")
    for f in ["Q12866.json", "Q60805.json", "summary.json", "residue_mapping.json"]:
        shutil.copy2(species / f, out / ("species_" + f))
    for f in ["P30530.json", "Q06418.json", "P07225.json", "Q14393.json", "tam_mapping.json"]:
        shutil.copy2(source / f, out / f)
    print(f"Prepared {len(entries)} structures and {len(patches)} exploratory patches in {out}")


if __name__ == "__main__":
    main()
