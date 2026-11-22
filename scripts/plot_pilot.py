#!/usr/bin/env python3
"""Plot native pilot geometry; these figures do not establish binding or agonism."""

import argparse
import json
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionException

from mertk.provenance import file_hash
from mertk.structural import rigid_fit, rmsd
from mertk.structures import read_design_chain

COLORS = {"bindcraft": "#007F86", "rfd4_proteina": "#BD5C27"}


def collect(campaign):
    reference, mapping = read_design_chain(campaign / "inputs/target.pdb", "A")
    records, skipped, seen = [], [], set()
    inputs = [("bindcraft", p) for p in sorted((campaign / "bindcraft/designs").rglob("*.pdb"))] + [
        ("rfd4_proteina", p)
        for p in sorted((campaign / "rfd4_proteina/inference_samples").rglob("*_generated.cif"))
    ]
    for arm, path in inputs:
        digest = file_hash(path)
        if digest in seen:
            skipped.append({"path": str(path.relative_to(campaign)), "reason": "identical file"})
            continue
        seen.add(digest)
        try:
            parser = MMCIFParser(QUIET=True) if path.suffix == ".cif" else PDBParser(QUIET=True)
            structure = parser.get_structure("pilot", str(path))
            if len(structure) != 1 or len(structure[0]) != 2:
                raise ValueError("Not a single-model, two-chain complex")
            chains = {c.id: read_design_chain(path, c.id)[0] for c in structure[0]}
            matches = [key for key, chain in chains.items() if chain.sequence == reference.sequence]
            if len(matches) != 1:
                raise ValueError(
                    "Cannot uniquely identify the complete experimental target sequence"
                )
            target_key = matches[0]
            binder_key = next(key for key in chains if key != target_key)
            target, binder = chains[target_key].backbone[:, 1], chains[binder_key].backbone[:, 1]
            rotation, translation = rigid_fit(target, reference.backbone[:, 1])
            distances = np.linalg.norm(target[:, None] - binder[None, :], axis=-1)
            records.append(
                {
                    "arm": arm,
                    "name": path.stem,
                    "path": str(path.relative_to(campaign)),
                    "sha256": digest,
                    "target_chain": target_key,
                    "binder_chain": binder_key,
                    "binder_sequence": chains[binder_key].sequence,
                    "target_ca_rmsd": rmsd(
                        target @ rotation + translation, reference.backbone[:, 1]
                    ),
                    "target_residues_with_ca_contact": int((distances <= 8).any(axis=1).sum()),
                    "binder_residues_with_ca_contact": int((distances <= 8).any(axis=0).sum()),
                    "ca_contact_pairs": int((distances <= 8).sum()),
                    "target": target @ rotation + translation,
                    "binder": binder @ rotation + translation,
                    "minimum_ca_distance": distances.min(axis=1),
                }
            )
        except (ValueError, KeyError, PDBConstructionException) as error:
            skipped.append({"path": str(path.relative_to(campaign)), "reason": str(error)})
    return records, skipped, mapping


def save(fig, directory, name):
    for extension in ("png", "svg"):
        fig.savefig(directory / f"{name}.{extension}", dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_af3(campaign, directory):
    scores = campaign / "assessment/af3/scores"
    if not (scores / "summary.json").exists():
        return False
    summary = json.loads((scores / "summary.json").read_text())
    if not summary.get("execution_verified"):
        return False
    rows = [json.loads(line) for line in (scores / "candidates.jsonl").read_text().splitlines()]
    valid = [row for row in rows if row["prediction_status"] == "complete"]
    if not valid:
        return False
    fig, ax = plt.subplots(figsize=(7, 5))
    for arm, color in COLORS.items():
        points = [row["metrics"] for row in valid if row["arm"] == arm]
        if points:
            ax.scatter(
                [p["pair_iptm"] for p in points],
                [p["binder_pose_rmsd"] for p in points],
                color=color,
                label=arm,
                s=70,
                alpha=0.8,
            )
    gates = summary["protocol"]["gates"]
    ax.axvline(gates["min_pair_iptm"], color="#7C8590", ls="--", lw=1)
    ax.axhline(gates["max_binder_pose_rmsd"], color="#7C8590", ls="--", lw=1)
    ax.set(
        xlim=(0, 1),
        xlabel="Shared AF3 target–binder pair ipTM",
        ylabel="Binder pose RMSD from its native design pose (Å)",
        title=f"Matched single-sequence AF3 pilot · {len(valid)}/{len(rows)} candidates scored\nStructural consistency only; no biological validation",
    )
    ax.legend(frameon=False)
    ax.grid(alpha=0.2)
    fig.tight_layout()
    save(fig, directory, "af3_confidence_pose")
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign", type=Path)
    args = parser.parse_args()
    directory = args.campaign / "figures"
    directory.mkdir(exist_ok=True)
    records, skipped, mapping = collect(args.campaign)
    summary = {
        "interpretation": "Native generator structures, including rejected/intermediate files; no shared AF3 assessment or biological validation. Counts are files, not independent successful designs.",
        "records": [
            {k: v for k, v in row.items() if k not in {"target", "binder", "minimum_ca_distance"}}
            for row in records
        ],
        "skipped": skipped,
    }
    (directory / "geometry.json").write_text(json.dumps(summary, indent=2) + "\n")
    if not records:
        (directory / "README.md").write_text(
            "No complete two-chain designs are available yet. No result figures were fabricated.\n"
        )
        print("No generated complexes to plot.")
        return
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    # Show up to three files from each arm; the JSON inventory includes every file.
    shown = [row for arm in COLORS for row in [r for r in records if r["arm"] == arm][:3]]
    fig = plt.figure(figsize=(12, 4 * ((len(shown) + 2) // 3)))
    for i, row in enumerate(shown):
        ax = fig.add_subplot((len(shown) + 2) // 3, 3, i + 1, projection="3d")
        for key, color, width in (("target", "#AEB7C2", 1), ("binder", COLORS[row["arm"]], 2)):
            points = row[key]
            ax.plot(*points.T, color=color, lw=width)
        points = np.concatenate([row["target"], row["binder"]])
        center, half = (points.min(0) + points.max(0)) / 2, np.ptp(points, axis=0).max() / 2 + 3
        ax.set(
            xlim=(center[0] - half, center[0] + half),
            ylim=(center[1] - half, center[1] + half),
            zlim=(center[2] - half, center[2] + half),
        )
        ax.set_box_aspect((1, 1, 1))
        ax.set_title(f"{row['arm']}\n{textwrap.fill(row['name'], 30)}", fontsize=8)
        ax.set_axis_off()
    fig.suptitle(
        "Native pilot complexes · gray MERTK, colored binder\nAligned to experimental 37NT; unvalidated predictions",
        fontsize=12,
    )
    fig.tight_layout()
    save(fig, directory, "complex_gallery")

    fig, ax = plt.subplots(figsize=(12, max(3, 0.5 * len(shown) + 2)))
    distances = np.stack([r["minimum_ca_distance"] for r in shown])
    image = ax.imshow(distances, aspect="auto", vmin=0, vmax=30, cmap="magma_r")
    ticks = np.linspace(0, len(mapping) - 1, 9, dtype=int)
    ax.set_xticks(ticks, [mapping[t]["author_number"] for t in ticks])
    ax.set_yticks(range(len(shown)), [f"{r['arm']}: {r['name']}" for r in shown], fontsize=8)
    ax.set(
        xlabel="Experimental MERTK author residue number",
        title="Proximity to the generated binder · Cα geometry only",
    )
    fig.colorbar(image, ax=ax, label="Nearest binder Cα distance (Å; capped at 30)")
    fig.tight_layout()
    save(fig, directory, "target_contact_profile")

    fig, ax = plt.subplots(figsize=(7, 5))
    for arm, color in COLORS.items():
        subset = [r for r in records if r["arm"] == arm]
        if subset:
            ax.scatter(
                [r["target_ca_rmsd"] for r in subset],
                [r["target_residues_with_ca_contact"] for r in subset],
                color=color,
                label=f"{arm} ({len(subset)} files)",
                s=60,
                alpha=0.8,
            )
    ax.set(
        xlabel="MERTK Cα RMSD from experimental 37NT after alignment (Å)",
        ylabel="Target residues within 8 Å of a binder Cα",
        title="Target preservation and geometric contact coverage\nNative files; no shared confidence or activity assessment",
    )
    ax.legend(frameon=False)
    ax.grid(alpha=0.2)
    fig.tight_layout()
    save(fig, directory, "geometry_overview")
    (directory / "README.md").write_text(
        "# Pilot geometry figures\n\n" + summary["interpretation"] + "\n\n"
        "The gallery/contact profile shows at most three files per arm, sorted by path; "
        "this is an inventory selection, not a ranking. The scatter shows all distinct "
        "readable complex files. Exact duplicate bytes are excluded. Relaxed structures "
        "and sequence variants can share one backbone/trajectory and are not independent replicates. "
        "`geometry.json` records file hashes, chain assignments, measured geometry and skipped files. "
        "Cα proximity is a coarse geometric proxy, not an affinity or agonism score.\n"
    )
    af3_plotted = plot_af3(args.campaign, directory)
    if af3_plotted:
        with (directory / "README.md").open("a") as handle:
            handle.write(
                "\nThe AF3 panel uses verified shared execution records; dashed lines show only "
                "the pair-confidence and pose gates (other structural gates are in the score report). "
                "This tiny, selected runtime cohort is not a generator-performance benchmark. "
                "AF3-derived summaries are subject to the [AF3 Output Terms of Use]"
                "(https://github.com/google-deepmind/alphafold3/blob/main/OUTPUT_TERMS_OF_USE.md).\n"
            )
    print(
        f"Wrote {3 + int(af3_plotted)} PNG/SVG figures from {len(records)} native complex files to {directory}"
    )


if __name__ == "__main__":
    main()
