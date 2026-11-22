#!/usr/bin/env python3
"""Reproduce the provisional Ig1 anchor dossier; AXL overlays are analogies only."""

import argparse
import json
from pathlib import Path

import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.PDB import MMCIFParser
from Bio.SeqUtils import seq1

from mertk.provenance import file_hash
from mertk.structural import rigid_fit, rmsd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("outputs/epitope_review_20261112"))
    args = parser.parse_args()
    native_path = Path("data/references/37NT/37NT.cif")
    axl_path = Path("data/references/pymol_review/2C5D.cif")
    sasa_path = Path("data/references/pymol_review/mertk_sasa.json")
    mapping_path = Path("outputs/species_review_20261111/residue_mapping.json")
    native = MMCIFParser(QUIET=True).get_structure("mertk", native_path)[0]
    axl = MMCIFParser(QUIET=True).get_structure("axl", axl_path)[0]
    sasa = json.loads(sasa_path.read_text())
    species = json.loads(mapping_path.read_text())
    epitopes = json.loads(Path("configs/course_epitopes.json").read_text())
    target = [r for r in native["A"] if r.id[0] == " " and r.id[1] <= 285]
    homolog = [r for r in axl["C"] if r.id[0] == " " and "CA" in r]
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score, aligner.extend_gap_score = -10, -0.5
    alignment = aligner.align(
        "".join(seq1(r.resname) for r in target),
        "".join(seq1(r.resname) for r in homolog),
    )[0]
    pairs = [
        (i, j)
        for (start, end), (other_start, other_end) in zip(*alignment.aligned)
        for i, j in zip(range(start, end), range(other_start, other_end))
        if 93 <= target[i].id[1] <= 186
    ]
    x = np.array([homolog[j]["CA"].coord for i, j in pairs])
    y = np.array([target[i]["CA"].coord for i, j in pairs])
    rotation, translation = rigid_fit(x, y)
    gas6 = {
        chain: np.array([a.coord for r in axl[chain] for a in r if a.element != "H"]) @ rotation
        + translation
        for chain in ["A", "B"]
    }
    glycans = np.array(
        [
            a.coord
            for r in native.get_residues()
            if r.resname in {"NAG", "MAN", "BMA"}
            for a in r
            if a.element != "H"
        ]
    )
    records = {}
    for key, epitope in epitopes.items():
        rows = []
        for n in epitope["hotspots"]:
            residue = native["A"][n]
            if not np.allclose(residue["CA"].coord, sasa[str(n)]["ca"], atol=0.001):
                raise ValueError("Solvent-area inventory does not match native coordinates")
            coordinates = np.array([a.coord for a in residue if a.element != "H"])
            row = {
                "author_number": n,
                "residue": residue.resname,
                "sasa_A2": sasa[str(n)]["sasa"],
                "observed_glycan_min_heavy_A": float(
                    np.linalg.norm(coordinates[:, None] - glycans[None, :], axis=-1).min()
                ),
                "species_alignment": next(r for r in species if r["human_position"] == n),
            }
            if key == "Ig1_side":
                homolog_residue = next(homolog[j] for i, j in pairs if target[i].id[1] == n)
                row["AXL_sequence_correspondence"] = {
                    "author_number": homolog_residue.id[1],
                    "residue": homolog_residue.resname,
                    "interpretation": "Homology correspondence, not measured selectivity",
                }
                row["AXL_analogy_GAS6_min_heavy_A_by_chain"] = {
                    chain: float(
                        np.linalg.norm(coordinates[:, None] - points[None, :], axis=-1).min()
                    )
                    for chain, points in gas6.items()
                }
            rows.append(row)
        records[key] = rows
    report = {
        "interpretation": "Provisional anchors; exposure and analogy do not establish binding, ligand compatibility, cooperativity or shedding protection",
        "source_sha256": {
            str(p): file_hash(p) for p in [native_path, axl_path, sasa_path, mapping_path]
        },
        "AXL_to_MERTK_Ig1_alignment": {
            "method": "Ig-pair sequence alignment BLOSUM62, gap open -10, extend -0.5; rigid fit Ig1 correspondences only",
            "correspondences": len(pairs),
            "ca_rmsd_A": rmsd(x @ rotation + translation, y),
            "residue_pairs_MERTK_AXL": [[target[i].id[1], homolog[j].id[1]] for i, j in pairs],
            "rotation": rotation.tolist(),
            "translation": translation.tolist(),
            "limitation": "Both GAS6 chains are inspected; this is not a MERTK–GAS6 or MERTK–PROS1 structure",
        },
        "epitopes": records,
    }
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "review.json").write_text(json.dumps(report, indent=2) + "\n")
    print(args.output / "review.json")


if __name__ == "__main__":
    main()
