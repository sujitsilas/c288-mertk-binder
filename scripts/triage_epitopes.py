#!/usr/bin/env python3
"""Rank candidate epitopes on the intact MERTK ectodomain before spending GPU time.

The November 12 campaigns chose both epitopes by eye, from solvent exposure, mouse
conservation and paralogue sequence differences. None of those predict *designability*,
and the resulting sites were poor on criteria that are cheap to measure first:

  - FNIIIa Patch A (I365/W367/A369) holds the receptor's best hydrophobic anchors, but
    its centroid sits 10.5 A from sequon N316 and every patch residue falls inside a
    glycan shadow. Anchor-level clearance to *resolved* sugars understated this, because
    four of the 14 sequons carry no modelled glycan in the insect-cell structure.
  - Ig1 K178/E183/V185 contains no F/Y/W/I/L/M residue at all, the anchor identities
    upstream guidance asks for. Its designs compensated with Met/Trp enrichment and a
    0.57 surface hydrophobicity in one trajectory, past the 0.35 filter; a second
    scored 0.32. These observations do not establish a causal explanation.

This screen scores every exposed-residue-centred patch on the intact A93-A484 reference
for anchor content, local concavity and whole-footprint clearance from all sequons, and
reports each site's separation from the residues implicated in natural-ligand binding so
a blocking and a non-blocking arm can be chosen deliberately.

Usage:
    .venv/bin/python scripts/triage_epitopes.py [--output outputs/epitope_triage.json]
"""

import argparse
import json
import math
import statistics
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser, ShrakeRupley

from mertk.provenance import file_hash

# Theoretical per-residue maxima in Gly-X-Gly (Tien et al., 2013, PLoS ONE 8:e80635).
MAX_SASA = {
    "ALA": 129,
    "ARG": 274,
    "ASN": 195,
    "ASP": 193,
    "CYS": 167,
    "GLN": 225,
    "GLU": 223,
    "GLY": 104,
    "HIS": 224,
    "ILE": 197,
    "LEU": 201,
    "LYS": 236,
    "MET": 224,
    "PHE": 240,
    "PRO": 159,
    "SER": 155,
    "THR": 172,
    "TRP": 285,
    "TYR": 263,
    "VAL": 174,
}
THREE_TO_ONE = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}
ANCHOR_RESIDUES = set("FYWILM")  # upstream guidance's preferred hotspot identities
SUGARS = {"NAG", "BMA", "MAN", "FUC", "GAL", "NDG", "A2G", "SIA"}
DOMAINS = [("Ig1", 93, 185), ("Ig2", 186, 285), ("FNIIIa", 286, 381), ("FNIIIb", 382, 484)]
# Human positions whose aligned mouse counterparts are functional ligand-dependent
# phagocytosis mutants; treated as a proximity flag, not an atomic ligand footprint.
LIGAND_IMPLICATED = [145, 147, 268, 274]

PATCH_RADIUS = 10.0  # approximate minibinder footprint radius
GLYCAN_EXCLUSION = 15.0  # conservative mammalian complex N-glycan reach from Asn
MIN_PATCH_RESIDUES = 6
EXPOSURE_CENTRE = 0.25  # relative SASA required of a patch-defining residue
EXPOSURE_MEMBER = 0.15
EXPOSURE_ANCHOR = 0.30


def domain_of(number):
    for name, low, high in DOMAINS:
        if low <= number <= high:
            return name
    return "outside"


def find_sequons(sequence_by_number):
    """N-X-S/T with X != P, over consecutively numbered observed residues."""
    numbers = sorted(sequence_by_number)
    found = []
    for index, number in enumerate(numbers[:-2]):
        if numbers[index + 1] != number + 1 or numbers[index + 2] != number + 2:
            continue
        asn, middle, hydroxyl = (sequence_by_number[number + offset] for offset in range(3))
        if asn == "N" and middle != "P" and hydroxyl in "ST":
            found.append(number)
    return found


def load(target_pdb, reference_cif):
    structure = PDBParser(QUIET=True).get_structure("mertk", str(target_pdb))
    model = structure[0]
    ShrakeRupley().compute(model, level="R")
    residues = [r for r in model["A"] if r.id[0] == " " and r.get_resname() in MAX_SASA]
    glycan_atoms = []
    if reference_cif and Path(reference_cif).exists():
        cif = MMCIFParser(QUIET=True).get_structure("cif", str(reference_cif))
        glycan_atoms = [
            a.coord
            for r in cif[0].get_residues()
            if r.get_resname().strip() in SUGARS
            for a in r
            if a.element != "H"
        ]
    return residues, glycan_atoms


def score_patches(residues, glycan_atoms):
    sequence = {r.id[1]: THREE_TO_ONE[r.get_resname()] for r in residues}
    sasa = {r.id[1]: r.sasa for r in residues}
    relative = {r.id[1]: r.sasa / MAX_SASA[r.get_resname()] for r in residues}
    alpha = {r.id[1]: r["CA"].coord for r in residues if "CA" in r}
    heavy = {r.id[1]: [a.coord for a in r if a.element != "H"] for r in residues}
    sequons = find_sequons(sequence)
    every_atom = np.array([c for atoms in heavy.values() for c in atoms])

    patches = []
    for centre in sorted(n for n in alpha if relative[n] > EXPOSURE_CENTRE):
        members = [
            n
            for n in alpha
            if relative[n] > EXPOSURE_MEMBER and math.dist(alpha[n], alpha[centre]) <= PATCH_RADIUS
        ]
        if len(members) < MIN_PATCH_RESIDUES:
            continue
        centroid = np.array(
            [statistics.mean(float(alpha[n][k]) for n in members) for k in range(3)]
        )
        anchors = [
            n for n in members if sequence[n] in ANCHOR_RESIDUES and relative[n] > EXPOSURE_ANCHOR
        ]
        shielded = sum(
            1
            for n in members
            if any(math.dist(alpha[n], alpha[s]) < GLYCAN_EXCLUSION for s in sequons)
        )
        nearest_sequon = min(
            ((float(np.linalg.norm(centroid - alpha[s])), s) for s in sequons),
            default=(99.0, None),
        )
        patches.append(
            {
                "centre": centre,
                "residue": sequence[centre],
                "domain": domain_of(centre),
                "members": members,
                "patch_sasa": round(sum(sasa[n] for n in members), 1),
                "anchors": [f"{sequence[n]}{n}" for n in anchors],
                "anchor_sasa": sum(sasa[n] for n in anchors),
                # Buried centroid means a groove; a flat beta-sheet face scores low.
                "neighbour_density": int(
                    np.count_nonzero(np.linalg.norm(every_atom - centroid, axis=1) < 12.0)
                ),
                "nearest_sequon_angstroms": round(nearest_sequon[0], 1),
                "nearest_sequon": nearest_sequon[1],
                "shielded_fraction": round(shielded / len(members), 2),
                "patch_min_observed_glycan_heavy_atom_A": (
                    round(
                        float(
                            np.linalg.norm(
                                np.array([a for n in members for a in heavy[n]])[:, None]
                                - np.asarray(glycan_atoms)[None, :],
                                axis=-1,
                            ).min()
                        ),
                        2,
                    )
                    if len(glycan_atoms)
                    else None
                ),
                "ligand_face_angstroms": round(
                    min(
                        min(math.dist(x, y) for x in heavy[centre] for y in heavy[implicated])
                        for implicated in LIGAND_IMPLICATED
                        if implicated in heavy
                    ),
                    1,
                ),
            }
        )

    densities = [p["neighbour_density"] for p in patches]
    low, high = min(densities), max(densities)
    best_anchor = max(p["anchor_sasa"] for p in patches) or 1.0
    for patch in patches:
        concavity = (patch["neighbour_density"] - low) / (high - low) if high > low else 0.0
        anchor = patch["anchor_sasa"] / best_anchor
        glycan = min(patch["nearest_sequon_angstroms"] / GLYCAN_EXCLUSION, 1.0) * (
            1 - patch["shielded_fraction"]
        )
        patch.update(
            anchor_score=round(anchor, 3),
            concavity_score=round(concavity, 3),
            glycan_score=round(glycan, 3),
            score=round(0.40 * anchor + 0.30 * concavity + 0.30 * glycan, 3),
        )
        patch.pop("anchor_sasa")
    patches.sort(key=lambda p: -p["score"])
    return patches, sequons, alpha


def deduplicate(patches, alpha, separation=12.0):
    kept = []
    for patch in patches:
        if any(
            math.dist(alpha[patch["centre"]], alpha[other["centre"]]) < separation for other in kept
        ):
            continue
        kept.append(patch)
    return kept


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--target-pdb", type=Path, default=Path("data/references/37NT/target.pdb"))
    parser.add_argument("--reference-cif", type=Path, default=Path("data/references/37NT/37NT.cif"))
    parser.add_argument("--output", type=Path, default=Path("outputs/epitope_triage.json"))
    parser.add_argument("--top", type=int, default=15)
    args = parser.parse_args()

    residues, glycan_atoms = load(args.target_pdb, args.reference_cif)
    patches, sequons, alpha = score_patches(residues, glycan_atoms)
    representatives = deduplicate(patches, alpha)

    print(f"Observed residues        : {len(residues)}")
    print(f"N-X-S/T sequons          : {len(sequons)} -> {sequons}")
    print(f"Resolved glycan atoms    : {len(glycan_atoms)}")
    print(f"Scored patches           : {len(patches)} ({len(representatives)} non-overlapping)")
    print()
    header = (
        f"{'rank':>4} {'centre':>7} {'domain':>7} {'score':>6} {'anchor':>7} "
        f"{'concav':>7} {'glycan':>7} {'sequon':>7} {'ligand':>7}  anchors"
    )
    print(header)
    print("-" * len(header))
    for rank, patch in enumerate(representatives[: args.top], start=1):
        print(
            f"{rank:>4} {patch['residue']}{patch['centre']:<6} {patch['domain']:>7} "
            f"{patch['score']:>6.3f} {patch['anchor_score']:>7.3f} "
            f"{patch['concavity_score']:>7.3f} {patch['glycan_score']:>7.3f} "
            f"{patch['nearest_sequon_angstroms']:>7.1f} "
            f"{patch['ligand_face_angstroms']:>7.1f}  "
            f"{','.join(patch['anchors'][:6]) or '-'}"
        )
    print()
    print("sequon  = centroid distance to nearest N-X-S/T attachment, A (higher is better)")
    print("ligand  = closest approach to a ligand-implicated residue, A; a low value")
    print("          flags proximity only; it does not establish ligand competition")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(
            {
                "target_pdb": str(args.target_pdb),
                "target_sha256": file_hash(args.target_pdb),
                "reference_cif_sha256": file_hash(args.reference_cif),
                "script_sha256": file_hash(__file__),
                "total_scored_patches": len(patches),
                "method": {
                    "patch_radius_angstroms": PATCH_RADIUS,
                    "glycan_exclusion_angstroms": GLYCAN_EXCLUSION,
                    "anchor_residues": sorted(ANCHOR_RESIDUES),
                    "ligand_implicated_residues": LIGAND_IMPLICATED,
                    "weights": {"anchor": 0.40, "concavity": 0.30, "glycan": 0.30},
                },
                "limitations": [
                    "Scores rank designability and accessibility hypotheses; they are not "
                    "binding, affinity or selectivity measurements.",
                    "Glycan exclusion uses a uniform radius from every sequon and does not model "
                    "mammalian glycoform branching, occupancy or heterogeneity.",
                    "Ligand proximity derives from aligned mouse functional mutants, not from a "
                    "solved human MERTK-GAS6 or MERTK-PROS1 complex.",
                    "Neighbour density is a heuristic, not a validated concavity or designability metric. "
                    "Overlapping patches are correlated; a centroid outside 15 A does not make a footprint glycan-clear.",
                ],
                "sequons": sequons,
                "patches": representatives,
                "all_patches": patches,
            },
            indent=2,
        )
        + "\n"
    )
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
