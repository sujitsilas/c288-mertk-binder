"""Plot day-adjusted MerTK-axis RNA effects from reconstructed sample-level counts."""

import argparse
import csv
import math
from pathlib import Path

from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[2]
FIGURES = ROOT / "assets/figures"
GENES = ("Mertk", "Gas6", "Pros1", "Axl", "Adam17", "Adam8", "Mmp8", "Mmp9", "Mmp12")
POPULATIONS = {"Mφ": "Macrophages", "Mono/MDM": "Mono/MDM"}
FIELDS = ("cell_type", "gene", "logFC", "CI.L", "CI.R", "adj.P.Val", "n_Burn", "n_Sham")
COLORS = {"Mφ": "#12675c", "Mono/MDM": "#76539c"}
INK = "#263b43"
SECONDARY = "#51676f"


def read_rows(path):
    selected = {}
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            key = row["cell_type"], row["gene"]
            if key[0] not in POPULATIONS or key[1] not in GENES:
                continue
            if key in selected:
                raise ValueError(f"Duplicate group and gene in {path}: {key}")
            selected[key] = {field: row[field] for field in FIELDS}
    expected = {(group, gene) for group in POPULATIONS for gene in GENES}
    if missing := expected - selected.keys():
        raise ValueError(f"Missing group and gene combinations in {path}: {sorted(missing)}")
    return selected


def draw(rows):
    figure = Figure(figsize=(8.4, 3.9), facecolor="white")
    FigureCanvasAgg(figure)
    axes = figure.subplots(ncols=2, sharex=True, sharey=True)
    figure.subplots_adjust(left=0.10, right=0.97, bottom=0.23, top=0.84, wspace=0.18)
    lower, upper = 0.0, 0.0
    y = list(range(len(GENES)))
    for ax, (group, title) in zip(axes, POPULATIONS.items(), strict=True):
        group_rows = [rows[group, gene] for gene in GENES]
        replicates = {(int(row["n_Burn"]), int(row["n_Sham"])) for row in group_rows}
        if len(replicates) != 1:
            raise ValueError(f"Inconsistent sample counts for {group}")
        ((burn, sham),) = replicates
        ax.set_title(
            f"{title}\n{burn} burn · {sham} sham samples", loc="left", color=INK, fontsize=10.5
        )
        ax.set_yticks(y, GENES, fontstyle="italic")
        ax.set_ylim(len(GENES) - 0.5, -0.5)
        ax.axvline(0, color="#8a9b9f", linewidth=1.0, linestyle="--", zorder=1)
        ax.axhline(3.5, color="#dce5e4", linewidth=1, zorder=1)
        color = COLORS[group]
        for index, row in enumerate(group_rows):
            effect, left, right, fdr = (float(row[key]) for key in FIELDS[2:6])
            if not all(math.isfinite(value) for value in (effect, left, right, fdr)):
                raise ValueError(f"Non-finite estimate in {row}")
            if not left <= effect <= right or not 0 <= fdr <= 1:
                raise ValueError(f"Invalid confidence interval or FDR in {row}")
            lower, upper = min(lower, left), max(upper, right)
            ax.errorbar(
                effect,
                index,
                xerr=[[effect - left], [right - effect]],
                marker="o",
                markersize=5.6,
                markerfacecolor=color if fdr < 0.05 else "white",
                markeredgecolor=color,
                ecolor=color,
                elinewidth=1.25,
                capsize=2.4,
                zorder=3,
            )
        ax.spines[["top", "left", "right"]].set_visible(False)
        ax.spines["bottom"].set_color("#a8b8ba")
        ax.tick_params(axis="y", length=0, labelsize=10, colors=INK, pad=8)
        ax.tick_params(axis="x", labelsize=9, colors=SECONDARY, length=3, pad=4)
        ax.set_axisbelow(True)
        ax.grid(axis="x", color="#e6edeb", linewidth=0.75)
    axes[0].set_xlim(math.floor(lower), math.ceil(upper))
    figure.supxlabel(
        "Burn − sham, day-adjusted log₂ fold change", color=INK, fontsize=10, y=0.12
    )
    figure.legend(
        handles=[
            Line2D(
                [],
                [],
                color=INK,
                marker="o",
                linestyle="none",
                markersize=5,
                label="Transcriptome FDR < 0.05",
            ),
            Line2D(
                [],
                [],
                color=INK,
                marker="o",
                markerfacecolor="white",
                linestyle="none",
                markersize=5,
                label="Transcriptome FDR ≥ 0.05",
            ),
            Line2D([], [], color=SECONDARY, linewidth=1.2, label="Pointwise 95% interval"),
        ],
        ncols=3,
        loc="lower center",
        bbox_to_anchor=(0.53, 0.018),
        frameon=False,
        fontsize=8.2,
    )
    return figure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=FIGURES / "burn_sham_mertk_axis.tsv")
    parser.add_argument("--output", type=Path, default=FIGURES / "burn_sham_mertk_axis.png")
    parser.add_argument(
        "--snapshot", type=Path, help="Save the plotted estimates as an aggregate source table"
    )
    args = parser.parse_args()
    rows = read_rows(args.input)
    figure = draw(rows)
    if args.snapshot:
        with args.snapshot.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows[group, gene] for group in POPULATIONS for gene in GENES)
    figure.savefig(args.output, dpi=240, facecolor="white")


if __name__ == "__main__":
    main()
