#!/usr/bin/env python3
"""Measure sensitivity of existing stalk alignments to a specified coordinate perturbation."""

import argparse
import json
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser

from mertk.cleavage_protector import STALK_RANGE
from mertk.provenance import file_hash
from mertk.structural import rigid_fit

CONFORMERS = ("base", "psi490_plus25", "psi490_minus25")
SEED = 2436
DRAWS = 200
COORDINATE_SD = 0.5
PROBE_RADIUS = 20.0
MEMBRANE_MARGIN = 5.0


def fit_sensitivity(points, noise):
    centroid = points.mean(axis=0)
    _, spread, axes = np.linalg.svd(points - centroid, full_matrices=False)
    probes = centroid + PROBE_RADIUS * axes[1:]
    offsets = []
    for shift in noise:
        rotation, translation = rigid_fit(points + shift, points)
        offsets.append(
            float(np.linalg.norm(probes @ rotation + translation - probes, axis=1).max())
        )
    quantiles = np.quantile(offsets, [0.5, 0.95, 1.0]).round(3)
    return {
        "coordinate_spread_A": spread.round(3).tolist(),
        "probe_displacement_A": dict(
            zip(("median", "p95", "max"), quantiles.tolist(), strict=True)
        ),
    }


def audit(references, screen_path):
    low, high = STALK_RANGE
    residues = range(low, high + 1)
    # Reuse the same noise draws so differences between conformers are paired.
    rng = np.random.default_rng(SEED)
    noise = rng.normal(0, COORDINATE_SD, size=(DRAWS, len(residues), 3))
    parser = PDBParser(QUIET=True)
    conformers, sources = {}, {}
    for name in CONFORMERS:
        path = references / f"{name}_full_context.pdb"
        chain = parser.get_structure(name, path)[0]["A"]
        points = np.array([chain[n]["CA"].coord for n in residues], dtype=float)
        conformers[name] = fit_sensitivity(points, noise)
        sources[path.name] = file_hash(path)

    screen = json.loads(screen_path.read_text())
    sources[screen_path.name] = file_hash(screen_path)
    fractions = [len(row["stalk_contacts"]) / max(1, len(row["contacts"])) for row in screen]
    clearance = np.array(
        [case["membrane_clearance_A"] for row in screen for case in row["conformers"]]
    )
    return {
        "settings": {
            "aligned_residues": STALK_RANGE,
            "draws": DRAWS,
            "seed": SEED,
            "per_coordinate_normal_sd_A": COORDINATE_SD,
            "probe_radius_A": PROBE_RADIUS,
            "summary_statistic": "Maximum displacement of two probes along the minor principal axes per draw.",
        },
        "conformers": conformers,
        "screen": {
            "poses": len(screen),
            "pass_all_three": sum(
                all(case["pass"] for case in row["conformers"]) for row in screen
            ),
            "median_fraction_stalk_contacts": round(float(np.median(fractions)), 3),
            "poses_with_fewer_than_half_contacts_on_stalk": sum(f < 0.5 for f in fractions),
            "poses_with_insufficient_or_off_stalk_contacts": sum(
                "Insufficient cleavage-site contact or predominantly off-stalk interface"
                in row["rejection_reasons"]
                for row in screen
            ),
            "pose_conformer_combinations": int(clearance.size),
            "below_membrane_margin": int((clearance < 0).sum()),
            "beyond_modeled_plane": int((clearance < -MEMBRANE_MARGIN).sum()),
            "minimum_clearance_relative_to_margin_A": round(float(clearance.min()), 3),
        },
        "sources_sha256": sources,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("references", type=Path)
    parser.add_argument("screen", type=Path)
    args = parser.parse_args()
    print(json.dumps(audit(args.references, args.screen), indent=2))


if __name__ == "__main__":
    main()
