#!/usr/bin/env python3
"""Build a portable PyMOL session using the Python shipped with PyMOL."""

import argparse
import json
from pathlib import Path

import numpy as np
import pymol2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    directory = args.directory.resolve()
    raw = directory / "structures"
    info = json.loads((directory / "review.json").read_text())
    with pymol2.PyMOL() as app:
        c = app.cmd
        c.viewport(1400, 900)
        c.bg_color("white")
        for setting, value in [
            ("orthoscopic", 1),
            ("antialias", 2),
            ("ray_shadows", 0),
            ("ambient", 0.45),
            ("specular", 0.15),
            ("cartoon_fancy_helices", 1),
            ("label_color", "black"),
            ("label_size", 14),
            ("float_labels", 0),
            ("label_outline_color", "white"),
            ("scene_buttons", 1),
            ("ray_opaque_background", 1),
            ("seq_view", 1),
        ]:
            c.set(setting, value)
        c.load(str(raw / "37NT.cif"), "MERTK_37NT")
        c.remove("MERTK_37NT and solvent")
        center = np.mean(c.get_coords("MERTK_37NT and polymer.protein"), axis=0)
        c.translate((-center).tolist(), "MERTK_37NT")
        c.hide("everything")
        c.show("cartoon", "MERTK_37NT")
        c.color("gray80", "MERTK_37NT")
        domains = [
            ("Ig1", 93, 186, "cyan"),
            ("Ig2", 197, 273, "marine"),
            ("FNIIIa", 286, 381, "teal"),
            ("FNIIIb", 386, 484, "salmon"),
        ]
        for name, lo, hi, color in domains:
            c.select(name, f"MERTK_37NT and polymer.protein and resi {lo}-{hi}")
            c.color(color, name)
        c.select("observed_glycans", "MERTK_37NT and resn NAG+MAN+BMA")
        c.show("sticks", "observed_glycans")
        c.color("gold", "observed_glycans")
        c.select(
            "potential_N_glycan_sites",
            "MERTK_37NT and polymer.protein and resi "
            + "+".join(map(str, info["potential_N_glycosylation_sequons"])),
        )
        c.select(
            "unmodeled_glycan_sites",
            "potential_N_glycan_sites and not resi "
            + "+".join(map(str, info["experimental_glycan_sites"])),
        )
        c.select(
            "human_mouse_changes",
            "MERTK_37NT and polymer.protein and resi "
            + "+".join(map(str, info["human_mouse_differences"])),
        )
        c.select("ligand_region_caution", "Ig1 or Ig2")
        c.select("crystal_contact_caution", "MERTK_37NT and polymer.protein and resi 444+449")
        for patch in info["patches"]:
            c.select(
                patch["id"],
                "MERTK_37NT and polymer.protein and resi " + "+".join(map(str, patch["anchors"])),
            )
        # These are separate molecules, never an invented MERTK-ligand complex.
        models = [
            ("1H30.cif", "GAS6_LG_experimental"),
            ("2C5D.cif", "AXL_GAS6_experimental"),
            ("1Z6C.cif", "PROS1_EGF34_experimental"),
            ("AF-P07225-F1-model_v6.cif", "PROS1_predicted"),
            ("AF-Q14393-F1-model_v6.cif", "GAS6_full_predicted"),
        ]
        for file, name in models:
            c.load(str(raw / file), name)
            c.remove(f"{name} and solvent")
            c.translate((-np.mean(c.get_coords(name, state=1), axis=0)).tolist(), name)
            c.hide("everything", name)
            c.show("cartoon", name)
        c.color("cyan", "GAS6_LG_experimental")
        c.color("lightblue", "AXL_GAS6_experimental and chain A+B")
        c.color("violet", "AXL_GAS6_experimental and chain C+D")
        for name, gla, egf, lg1, lg2 in [
            ("PROS1_predicted", (42, 87), (117, 283), (299, 475), (484, 666)),
            ("GAS6_full_predicted", (53, 94), (116, 278), (298, 470), (477, 670)),
        ]:
            c.color("gray70", name)
            for (lo, hi), color in [
                (gla, "green"),
                (egf, "marine"),
                (lg1, "yelloworange"),
                (lg2, "orange"),
            ]:
                c.color(color, f"{name} and resi {lo}-{hi}")
            c.hide("cartoon", f"{name} and resi 1-{gla[0] - 1}")
        c.color("marine", "PROS1_EGF34_experimental")
        for i, d in enumerate(info["rfd4_diagnostics"]):
            name = f"RFD4_rep{i}_failed_geometry"
            c.load(str(raw / d["file"]), name)
            xyz = (
                c.get_coords(name) @ np.asarray(d["rotation"])
                + np.asarray(d["translation"])
                - center
            )
            c.load_coords(xyz, name)
            c.hide("everything", name)
            c.show("cartoon", name)
            c.color("magenta" if i == 0 else "orange", f"{name} and chain A")
            c.color("lime", f"{name} and chain B")

        def scene(name, objects, message, selection=None, image=False, orient=True):
            c.disable("all")
            for obj in objects:
                c.enable(obj)
            c.set("grid_mode", 0)
            c.deselect()
            if orient:
                c.orient(selection or "enabled")
            else:
                c.reset()
            c.zoom(selection or "enabled", buffer=9)
            c.clip("slab", 300)
            c.refresh()
            c.scene(name, "store", message=message)
            if image:
                c.png(
                    str(directory / (name + ".png")),
                    width=1800,
                    height=1100,
                    dpi=150,
                    ray=1,
                    quiet=1,
                )

        scene(
            "01_MERTK_domains",
            ["MERTK_37NT"],
            "Experimental 37NT: Ig1 cyan, Ig2 blue, FNIIIa teal, FNIIIb salmon; observed glycans gold. Stalk and membrane absent.",
            image=True,
        )
        overview = []
        for x, source, name in [
            (-110, "MERTK_37NT", "view_MERTK"),
            (0, "GAS6_LG_experimental", "view_GAS6_LG"),
            (110, "PROS1_predicted", "view_PROS1_predicted"),
        ]:
            c.create(name, source, source_state=1, target_state=1)
            xyz = c.get_coords(name)
            xyz = xyz - xyz.mean(0)
            _, _, axes = np.linalg.svd(xyz, full_matrices=False)
            # Longest principal axis vertical; each molecule gets an independent orientation.
            rot = axes[[1, 0, 2]].T
            if np.linalg.det(rot) < 0:
                rot[:, 2] *= -1
            xyz = xyz @ rot + [x, 0, 0]
            c.load_coords(xyz, name)
            overview.append(name)
        scene(
            "02_three_molecules",
            overview,
            "Separate molecules, left to right: experimental MERTK, experimental GAS6 LG domains, predicted PROS1. Positions and orientations are for display only; no docking pose is implied.",
            image=True,
            orient=False,
        )
        c.show("surface", "MERTK_37NT and polymer.protein")
        c.set("transparency", 0.75, "MERTK_37NT")
        c.show("spheres", "unmodeled_glycan_sites and name ND2")
        c.set("sphere_scale", 0.55, "unmodeled_glycan_sites")
        c.color("yellow", "unmodeled_glycan_sites")
        for patch in info["patches"]:
            c.hide("labels", "all")
            c.hide("sticks", "MERTK_37NT and polymer.protein")
            c.color("gray80", "MERTK_37NT and polymer.protein")
            c.color("magenta", patch["id"])
            c.show("sticks", patch["id"])
            for n, label in zip(patch["anchors"], patch["anchor_labels"]):
                c.label(f"MERTK_37NT and polymer.protein and resi {n} and name CA", repr(label))
            scene(
                "03_patch_A" if patch["id"].startswith("A_") else "04_patch_B",
                ["MERTK_37NT"],
                "PROPOSED anchors only. " + patch["rationale"] + " " + patch["limitations"],
                selection="MERTK_37NT and (" + patch["id"] + " around 12)",
                image=True,
            )
        c.hide("labels", "all")
        c.hide("surface", "all")
        c.hide("spheres", "all")
        c.color("gray80", "MERTK_37NT and polymer.protein")
        c.color("magenta", "human_mouse_changes")
        c.show("surface", "MERTK_37NT and polymer.protein")
        c.set("transparency", 0.1, "MERTK_37NT")
        scene(
            "05_mouse_changes",
            ["MERTK_37NT"],
            "Magenta residues differ from canonical mouse MERTK. Conservation does not establish cross-reactive binding.",
            image=True,
        )
        scene(
            "06_AXL_GAS6",
            ["AXL_GAS6_experimental"],
            "Experimental 2C5D: GAS6 light blue, AXL violet. This is AXL, NOT MERTK; do not transfer the pose or contacts as established MERTK data.",
            image=True,
        )
        scene(
            "07_PROS1_domains",
            ["PROS1_predicted"],
            "Predicted PROS1: Gla green, EGF blue, LG1 yellow-orange, LG2 orange. Mature sequence shown; no glycans, gamma-carboxylation or validated receptor pose.",
            image=True,
        )
        c.spectrum("b", "red_yellow_green_cyan_blue", "PROS1_predicted", minimum=0, maximum=100)
        scene(
            "08_PROS1_confidence",
            ["PROS1_predicted"],
            "Predicted PROS1 colored by pLDDT (red low, blue high); domain arrangement is not experimentally established.",
        )
        scene(
            "09_GAS6_full",
            ["GAS6_full_predicted"],
            "Predicted full GAS6: Gla green, EGF blue, LG1 yellow-orange, LG2 orange. Experimental GAS6 LG domains are available in a separate object.",
        )
        scene(
            "10_PROS1_EGF_NMR",
            ["PROS1_EGF34_experimental"],
            "Experimental 1Z6C, EGF3-4 fragment only. State 1 of a 20-model NMR ensemble; receptor-binding LG domains are absent.",
        )
        c.hide("surface", "all")
        c.hide("labels", "all")
        c.color("gray70", "MERTK_37NT and polymer.protein")
        for i in range(len(info["rfd4_diagnostics"])):
            scene(
                f"11_RFD4_rep{i}",
                ["MERTK_37NT", f"RFD4_rep{i}_failed_geometry"],
                "Experimental MERTK gray; generated target colored; binder green. Global sequence-matched target alignment. Failed geometry, not a validated binder.",
                image=True,
            )
        c.scene("01_MERTK_domains", "recall", animate=0)
        c.save(str(directory / "MERTK_ligands_hotspots.pse"))
        (directory / "pymol_build.json").write_text(
            json.dumps(
                {
                    "pymol_version": c.get_version(),
                    "scenes": c.get_scene_list(),
                    "objects": c.get_names("objects"),
                },
                indent=2,
            )
            + "\n"
        )
        print("Saved session with", len(c.get_scene_list()), "scenes")


if __name__ == "__main__":
    main()
