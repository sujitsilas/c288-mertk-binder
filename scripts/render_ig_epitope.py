#!/usr/bin/env python3
"""Render the native experimental receptor and new anchors with PyMOL's Python."""

import argparse
import shutil
from pathlib import Path

import numpy as np
import pymol2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("outputs/epitope_review_20261112"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    source = Path("data/references/37NT/37NT.cif")
    shutil.copy2(source, args.output / "37NT.cif")
    with pymol2.PyMOL() as app:
        c = app.cmd
        c.viewport(1000, 1000)
        c.load(str(source.resolve()), "MERTK_37NT")
        c.remove("solvent")
        c.bg_color("white")
        for setting, value in [
            ("orthoscopic", 1),
            ("antialias", 2),
            ("ray_shadows", 0),
            ("label_color", "black"),
            ("label_size", 24),
            ("label_outline_color", -1),
            ("float_labels", 0),
            ("two_sided_lighting", 1),
            ("scene_buttons", 1),
        ]:
            c.set(setting, value)
        c.hide("everything")
        c.show("cartoon", "polymer.protein")
        c.color("gray80", "polymer.protein")
        for name, extent, color in [
            ("Ig1", "93-186", "cyan"),
            ("Ig2", "197-273", "marine"),
            ("FNIIIa", "286-381", "teal"),
            ("FNIIIb", "386-484", "salmon"),
        ]:
            c.select(name, f"MERTK_37NT and polymer.protein and resi {extent}")
            c.color(color, name)
        c.select("observed_glycans", "MERTK_37NT and resn NAG+MAN+BMA")
        c.show("sticks", "observed_glycans")
        c.color("gold", "observed_glycans")
        for name, residues, color in [
            ("Ig1_anchors", "178+183+185", "violet"),
            ("Patch_A", "365+367+369", "magenta"),
        ]:
            c.select(name, f"MERTK_37NT and polymer.protein and resi {residues}")
            c.show("sticks", name)
            c.color(color, name)
        c.orient("MERTK_37NT")
        c.zoom("MERTK_37NT", 10, complete=1)
        c.label("Ig1_anchors and resi 183 and name CA", '"Ig1 anchors"')
        c.label("Patch_A and resi 367 and name CA", '"Patch A"')
        c.deselect()
        c.scene("01_two_exploratory_epitopes", "store")
        c.png(str(args.output / "two_epitopes.png"), width=1000, height=1000, dpi=180, ray=1)
        c.show("surface", "Ig1")
        c.set("transparency", 0.2, "Ig1")
        c.viewport(1100, 900)
        c.label("all", '""')
        anchor = c.get_coords("Ig1_anchors").mean(0)
        center = c.get_coords("Ig1").mean(0)
        z = anchor - center
        z /= np.linalg.norm(z)
        ref = np.array([0.0, 0.0, 1.0])
        if abs(z @ ref) > 0.9:
            ref = np.array([0.0, 1.0, 0.0])
        x = np.cross(ref, z)
        x /= np.linalg.norm(x)
        y = np.cross(z, x)
        rotation = np.array([x, y, z])
        for n, label, offset in [(178, "K178", -9), (183, "E183", 0), (185, "V185", 9)]:
            point = anchor + z * 15 + y * 15 + x * offset
            c.pseudoatom(f"label_{n}", pos=point.tolist(), label=label)
            c.hide("nonbonded", f"label_{n}")
            c.color("black", f"label_{n}")
            c.distance(
                f"guide_{n}",
                f"label_{n}",
                f"MERTK_37NT and polymer.protein and resi {n} and name CA",
            )
            c.hide("labels", f"guide_{n}")
            c.set("dash_color", "gray50", f"guide_{n}")
            c.set("dash_radius", 0.025, f"guide_{n}")
        c.set_view([*rotation.T.flatten(), 0, 0, -100, *anchor, 50, 150, 0])
        c.zoom("Ig1 or label_*", 9, complete=1)
        c.clip("slab", 200)
        c.deselect()
        c.scene("02_Ig1_side_surface_hypothesis", "store")
        c.png(str(args.output / "Ig1_anchors.png"), width=1100, height=900, dpi=180, ray=1)
        c.save(str(args.output / "MERTK_Ig1_review.pse"))
    (args.output / "README.md").write_text(
        "# Exploratory Ig1 binder surface\n\n"
        "Open `MERTK_Ig1_review.pse` in PyMOL. Scene 01 shows both epitope hypotheses; "
        "scene 02 shows Ig1 K178/E183/V185. Gold sticks are observed glycans. "
        "Only the experimental unliganded receptor is displayed. Anchors are not validated hotspots. "
        "`review.json` records the separate AXL-based ligand analogy and its limitations. "
        "These images do not establish ligand compatibility, agonism, cooperativity or protection "
        "from shedding. The full stalk and mammalian glycan envelopes are unresolved.\n"
    )
    print(args.output / "MERTK_Ig1_review.pse")


if __name__ == "__main__":
    main()
