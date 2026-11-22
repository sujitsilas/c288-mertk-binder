"""Evidence and ranking for exploratory human MerTK stalk poses."""

import csv
import re

import numpy as np


def disable_score_cutoffs(settings):
    """Preserve BindCraft's metric schema, including per-amino-acid settings."""
    for item in settings.values():
        if "threshold" in item:
            item["threshold"] = None
        else:
            disable_score_cutoffs(item)
    return settings


def bond_distances(binder_xyz, receptor):
    return [
        float(np.linalg.norm(binder_xyz - receptor[n][atom].coord, axis=1).min())
        for n, atom in [(490, "C"), (491, "N")]
    ]


def rank_key(row):
    # Distances from an early-stopped trajectory can be artificially short in a clash.
    tiers = {"Sequence refold": 0, "Scored backbone": 1}
    return (
        tiers.get(row["stage"], 2),
        max(row["cleavage_atom_min_distances_A"]),
        -row["iptm"] if row["iptm"] is not None else float("inf"),
        row["name"],
    )


def saved_poses(designs):
    """Yield each saved candidate once, excluding duplicate relaxed backbone copies."""
    for folder, table, stage, score_prefix in [
        ("Accepted", "final_design_stats.csv", "Sequence refold", "Average_"),
        ("Trajectory", "trajectory_stats.csv", "Scored backbone", ""),
        ("Trajectory/Clashing", None, "Clash flagged", ""),
        ("Trajectory/LowConfidence", None, "Low confidence", ""),
    ]:
        stats = {}
        if table and (designs / table).is_file():
            with (designs / table).open() as stream:
                stats = {r["Design"]: r for r in csv.DictReader(stream)}
        for path in sorted((designs / folder).glob("*.pdb")):
            name = re.sub(r"_model\d+$", "", path.stem) if folder == "Accepted" else path.stem
            row = stats.get(name, {})
            # An independently scored sequence must be traceable to its final statistics.
            if folder == "Accepted" and not row:
                yield path, {"error": f"No final statistics for {path.name}"}
                continue
            yield (
                path,
                dict(
                    name=name,
                    stage=stage if row or not table else "Unscored backbone",
                    expected_sequence=row.get("Sequence"),
                    iptm=float(row[score_prefix + "i_pTM"])
                    if row.get(score_prefix + "i_pTM")
                    else None,
                    plddt=100 * float(row[score_prefix + "pLDDT"])
                    if row.get(score_prefix + "pLDDT")
                    else None,
                ),
            )
