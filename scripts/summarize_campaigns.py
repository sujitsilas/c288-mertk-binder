#!/usr/bin/env python3
"""Create auditable metrics and a static diagnostic figure from native campaign files."""

import argparse
import csv
import json
from pathlib import Path

from mertk.campaign_metrics import summarize


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaigns", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, default=Path("outputs/analysis/20261112"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    records = [summarize(p) for p in args.campaigns]
    (args.output / "campaign_metrics.json").write_text(json.dumps(records, indent=2) + "\n")
    fields = [
        "campaign_id",
        "attempts_started",
        "attempts_completed",
        "attempts_censored",
        "trajectories_scored",
        "saved_sequence_proposals",
        "base_af2_rejected_variants",
        "accepted_variants",
        "accepted_trajectories",
    ]
    with (args.output / "campaign_metrics.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows({k: row[k] for k in fields} for row in records)
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.4), layout="constrained")
    colours = ["#276b94", "#cb6944", "#4d967b"]
    names = [r.get("epitope") or r["campaign_id"] for r in records]
    categories = ["Started", "Scored\nbackbones", "Accepted\nbackbones"]
    x = np.arange(3)
    width = 0.8 / len(records)
    for i, row in enumerate(records):
        values = [row["attempts_started"], row["trajectories_scored"], row["accepted_trajectories"]]
        bars = axes[0].bar(
            x - 0.4 + width / 2 + i * width,
            values,
            width,
            color=colours[i % len(colours)],
            label=names[i],
        )
        axes[0].bar_label(bars, padding=3, fontsize=9)
    axes[0].set_xticks(x, categories)
    axes[0].set_ylabel("Independent trajectories")
    axes[0].set_title("Generation funnel")
    axes[0].legend(frameon=False, fontsize=9)
    failure_keys = ["pLDDT", "i_pTM", "i_pAE"]
    for i, row in enumerate(records):
        values = [row["failure_counts_overlapping"].get(k, 0) for k in failure_keys]
        bars = axes[1].bar(
            x - 0.4 + width / 2 + i * width, values, width, color=colours[i % len(colours)]
        )
        axes[1].bar_label(bars, padding=3, fontsize=9)
    axes[1].set_xticks(x, ["pLDDT", "interface\npTM", "interface\nPAE"])
    axes[1].set_ylabel("Recorded failures (overlap allowed)")
    axes[1].set_title("Sequence-validation bottleneck")
    for i, row in enumerate(records):
        for metric in row["trajectory_metrics"]:
            axes[2].scatter(
                float(metric["i_pTM"]),
                float(metric["Surface_Hydrophobicity"]),
                color=colours[i % len(colours)],
                s=70,
                edgecolors="white",
            )
    axes[2].axhline(0.35, linestyle="--", color="#666666", linewidth=1)
    axes[2].set_xlabel("Hallucinated trajectory interface pTM")
    axes[2].set_ylabel("Trajectory surface hydrophobicity")
    axes[2].set_title("Backbone-stage diagnostics")
    axes[2].text(
        0.03,
        0.95,
        "Dashed: later default filter = 0.35\nThese are not accepted sequence metrics",
        transform=axes[2].transAxes,
        va="top",
        fontsize=8,
    )
    axes[2].set_xlim(0, 1)
    axes[2].set_ylim(0, 0.7)
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", alpha=0.15)
        ax.set_axisbelow(True)
    fig.suptitle(
        "MERTK pilot: no accepted designs; confidence before sequence redesign is insufficient",
        fontsize=12,
    )
    fig.savefig(args.output / "pilot_metrics.png", dpi=180)
    fig.savefig(args.output / "pilot_metrics.svg")
    plt.close(fig)
    print(json.dumps([{k: r[k] for k in fields} for r in records], indent=2))


if __name__ == "__main__":
    main()
