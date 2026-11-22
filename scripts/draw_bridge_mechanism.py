#!/usr/bin/env python3
"""Draw a conceptual bridge/shedding map and a plate of experimental partial complexes."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle

ROOT = Path(__file__).resolve().parents[1]
COLORS = {
    "PS": "#A94F87",
    "Gla": "#31844B",
    "EGF": "#447CB1",
    "LG": "#C67D27",
    "Annexin": "#268F95",
    "MERTK": "#336A91",
    "muted": "#526372",
    "dark": "#21323D",
}


def block(ax, x, y, width, label, color, height=0.15):
    ax.add_patch(Rectangle((x, y), width, height, facecolor=color, edgecolor="white", linewidth=1))
    ax.text(
        x + width / 2,
        y + height / 2,
        label,
        ha="center",
        va="center",
        color="white",
        fontsize=11,
        weight="bold",
    )


def arrow(ax, x1, y1, x2, y2, **kwargs):
    ax.annotate(
        "",
        xy=(x2, y2),
        xytext=(x1, y1),
        arrowprops={"arrowstyle": "->", "lw": 1.5, "color": COLORS["muted"], **kwargs},
    )


def clean(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--bundle", type=Path, default=ROOT / "outputs/pymol/bridge_review_20261112"
    )
    args = parser.parse_args()
    out = ROOT / "assets/figures"
    out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {"font.family": "DejaVu Sans", "text.color": COLORS["dark"], "font.size": 11}
    )
    fig, axes = plt.subplots(3, 1, figsize=(12, 10.3), gridspec_kw={"height_ratios": [1, 1.1, 1]})
    fig.subplots_adjust(left=0.045, right=0.97, top=0.92, bottom=0.06, hspace=0.25)
    fig.suptitle(
        "PS recognition, receptor engagement and receptor preservation",
        fontsize=18,
        weight="bold",
        x=0.045,
        ha="left",
        y=0.975,
    )
    for ax in axes:
        clean(ax)
    ax = axes[0]
    ax.text(
        0,
        0.94,
        "A  Natural bridge: either GAS6 or PROS1 supplies both functions",
        fontsize=13,
        weight="bold",
    )
    for y in [0.225, 0.365, 0.505, 0.645]:
        ax.add_patch(Circle((0.055, y), 0.014, color=COLORS["PS"]))
        ax.plot([0.005, 0.04], [y - 0.025, y], color="#AAA9A2", lw=2)
        ax.plot([0.005, 0.04], [y + 0.025, y], color="#AAA9A2", lw=2)
    ax.text(0.005, 0.14, "Apoptotic-cell\nmembrane: PS exposed", fontsize=10, ha="left")
    ax.text(0.095, 0.76, "Ca²⁺", color=COLORS["Gla"], weight="bold")
    block(ax, 0.15, 0.43, 0.09, "Gla", COLORS["Gla"])
    block(ax, 0.265, 0.43, 0.13, "EGF spacer", COLORS["EGF"])
    block(ax, 0.42, 0.43, 0.15, "LG1 / LG2", COLORS["LG"])
    for left, right in [(0.075, 0.15), (0.24, 0.265), (0.395, 0.42), (0.57, 0.635)]:
        ax.plot([left, right], [0.505, 0.505], color=COLORS["muted"], lw=1.5)
    block(ax, 0.635, 0.43, 0.12, "Ig1 / Ig2", COLORS["MERTK"])
    block(ax, 0.765, 0.43, 0.12, "FNIIIa / b", COLORS["MERTK"])
    ax.plot([0.89, 0.91], [0.505, 0.505], color=COLORS["MERTK"], lw=2)
    ax.plot([0.91, 0.91], [0.22, 0.8], color="#AAA9A2", lw=6)
    ax.plot([0.91, 0.965], [0.505, 0.505], color=COLORS["MERTK"], lw=2)
    ax.text(0.69, 0.73, "MERTK on macrophage", fontsize=11)
    ax.text(0.15, 0.25, "Membrane binding", color=COLORS["Gla"], fontsize=10)
    ax.text(0.42, 0.25, "Receptor binding", color=COLORS["LG"], fontsize=10)
    ax.text(
        0.65,
        0.04,
        "Productive organization → kinase signaling → engulfment",
        fontsize=10,
        ha="center",
    )

    ax = axes[1]
    ax.text(
        0,
        0.94,
        "B  Swap the membrane-binding module; retain a receptor-binding module",
        fontsize=13,
        weight="bold",
    )
    block(ax, 0.01, 0.62, 0.2, "Annexin A5", COLORS["Annexin"])
    block(ax, 0.235, 0.62, 0.1, "Gla", COLORS["Gla"])
    ax.plot([0.21, 0.235], [0.695, 0.695], color=COLORS["muted"], lw=1.5)
    ax.text(
        0.38,
        0.695,
        "PS binding + PS binding; known MERTK-binding module absent",
        fontsize=11,
        va="center",
    )
    block(ax, 0.01, 0.29, 0.2, "Annexin A5", COLORS["Annexin"])
    block(ax, 0.265, 0.29, 0.22, "PROS1 LG1 / LG2", COLORS["LG"])
    ax.plot([0.21, 0.265], [0.365, 0.365], color=COLORS["muted"], lw=1.5, ls="--")
    ax.text(
        0.55,
        0.365,
        "Bridge hypothesis: PS recognition + receptor engagement",
        fontsize=11,
        va="center",
    )
    ax.text(
        0.01,
        0.04,
        "Domain boundaries, spacing and orientation require testing; PS binding alone does not establish agonism.",
        fontsize=10,
        color=COLORS["muted"],
    )

    ax = axes[2]
    ax.text(
        0,
        0.94,
        "C  Prevent shedding ≠ reattach an ectodomain that has already been released",
        fontsize=13,
        weight="bold",
    )
    block(ax, 0.01, 0.51, 0.20, "Intact MERTK", COLORS["MERTK"])
    arrow(ax, 0.225, 0.585, 0.37, 0.585)
    ax.text(0.295, 0.73, "ADAM17", ha="center", fontsize=10)
    block(ax, 0.385, 0.51, 0.235, "Soluble ectodomain", "#71889A")
    ax.text(0.645, 0.585, "+ membrane remnant", va="center", fontsize=11)
    ax.text(0.01, 0.23, "Before cleavage", weight="bold", fontsize=11)
    ax.text(
        0.01,
        0.08,
        "Test receptor-specific protection\nand useful activation separately.",
        fontsize=10,
    )
    ax.text(0.385, 0.23, "After cleavage", weight="bold", fontsize=11)
    ax.text(
        0.385,
        0.08,
        "No ordinary extracellular binder reconnects the receptor.\nRecovery needs intact receptor return or replacement.",
        fontsize=10,
    )
    fig.text(
        0.045,
        0.015,
        "Conceptual map, not atomic docking, active-complex stoichiometry or measured membrane spacing. Sources: Lew 2014; Thorp 2011; Cai 2016.",
        fontsize=9,
        color=COLORS["muted"],
    )
    fig.savefig(out / "bridge_first_principles.png", dpi=180, facecolor="white")
    fig.savefig(out / "bridge_first_principles.svg", facecolor="white")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
    fig.subplots_adjust(left=0.02, right=0.98, top=0.82, bottom=0.21, wspace=0.04)
    for ax, filename, title in zip(
        axes,
        ["04_Gla_lysoPS.png", "07_Annexin_contact_detail.png"],
        ["A  Gla–Ca–lysoPS: experimental surrogate", "B  Annexin V–Ca–headgroup: experimental"],
        strict=True,
    ):
        ax.imshow(mpimg.imread(args.bundle / filename))
        ax.axis("off")
        ax.set_title(title, fontsize=12, loc="left", pad=12)
    fig.suptitle(
        "Two ways to recognize an anionic lipid membrane",
        x=0.025,
        ha="left",
        y=0.97,
        fontsize=17,
        weight="bold",
    )
    fig.text(
        0.025,
        0.14,
        "1NL2: bovine prothrombin Gla + lysoPS analogue.\nCGU oxygens red; Ca green; lipid carbons magenta.\nNot a GAS6/PROS1 structure.",
        fontsize=10,
        linespacing=1.5,
    )
    fig.text(
        0.52,
        0.14,
        "1A8A: rat Annexin V + glycerophosphoserine.\nCa green; headgroup carbons magenta; no acyl tails.\nDashes: deposited Ca–O distances ≤3.2 Å.",
        fontsize=10,
        linespacing=1.5,
    )
    fig.text(
        0.025,
        0.015,
        "Our renderings of deposited coordinates; no invented ternary complex, fusion pose or neutrophil bilayer. See docs/first_principles_bridge.md for sources.",
        fontsize=9,
        color=COLORS["muted"],
    )
    fig.savefig(out / "lipid_binding_structural_evidence.png", dpi=180, facecolor="white")
    plt.close(fig)
    print("Saved mechanism map and experimental structural plate to", out)


if __name__ == "__main__":
    main()
