#!/usr/bin/env python3
"""Build first-principles PyMOL scenes; never fabricate a bound MERTK bridge."""

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import pymol2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    out = args.directory.resolve()
    raw = out / "structures"
    with pymol2.PyMOL() as app:
        c = app.cmd
        c.viewport(1600, 1000)
        c.bg_color("white")
        for setting, value in [
            ("orthoscopic", 1),
            ("antialias", 2),
            ("ray_shadows", 0),
            ("ambient", 0.45),
            ("specular", 0.15),
            ("label_color", "black"),
            ("label_size", 17),
            ("label_outline_color", "white"),
            ("scene_buttons", 1),
            ("stick_radius", 0.17),
        ]:
            c.set(setting, value)
        objects = [
            ("37NT.cif", "MERTK_experimental"),
            ("AF-Q14393-F1-model_v6.cif", "GAS6_predicted"),
            ("AF-P07225-F1-model_v6.cif", "PROS1_predicted"),
            ("1NL2.cif", "bovine_Gla_lysoPS_experimental"),
            ("1A8A.cif", "rat_Annexin_headgroup_experimental"),
            ("1ANX.cif", "human_Annexin_experimental"),
            ("2C5D.cif", "AXL_GAS6_experimental"),
            ("P5S_ideal.sdf", "PS_ideal_chemical_model"),
        ]
        for file, name in objects:
            c.load(str(raw / file), name)
            c.remove(f"{name} and solvent")
            c.hide("everything", name)
            c.show("cartoon", f"{name} and polymer.protein")
            c.color("gray75", name)
        c.remove("human_Annexin_experimental and not chain A")
        c.alter("PS_ideal_chemical_model", "resn='P5S'")
        c.remove("PS_ideal_chemical_model and hydro")
        c.remove("GAS6_predicted and resi 1-30")
        c.remove("PROS1_predicted and resi 1-41")
        for name, gla, egf, lg1, lg2 in [
            ("GAS6_predicted", (53, 94), (116, 278), (298, 470), (477, 670)),
            ("PROS1_predicted", (42, 87), (117, 283), (299, 475), (484, 666)),
        ]:
            for (lo, hi), color in [
                (gla, "forest"),
                (egf, "marine"),
                (lg1, "orange"),
                (lg2, "yelloworange"),
            ]:
                c.color(color, f"{name} and resi {lo}-{hi}")
        for lo, hi, color in [
            (93, 186, "cyan"),
            (197, 273, "marine"),
            (286, 381, "teal"),
            (386, 484, "salmon"),
        ]:
            c.color(color, f"MERTK_experimental and polymer.protein and resi {lo}-{hi}")
        c.show("sticks", "MERTK_experimental and resn NAG+MAN+BMA")
        c.color("gold", "MERTK_experimental and resn NAG+MAN+BMA")
        for name in ["rat_Annexin_headgroup_experimental", "human_Annexin_experimental"]:
            c.color("teal", f"{name} and polymer.protein")
            c.show("spheres", f"{name} and elem Ca")
        g = "bovine_Gla_lysoPS_experimental"
        c.hide("cartoon", g)
        c.show("cartoon", f"{g} and polymer.protein and resi 1-46")
        c.color("forest", f"{g} and polymer.protein")
        c.color("red", f"{g} and resn CGU and elem O")
        c.show("sticks", f"{g} and resn CGU+LPS")
        c.show("spheres", f"{g} and elem Ca")
        c.show("sticks", "rat_Annexin_headgroup_experimental and resn GSE")
        c.show("sticks", "PS_ideal_chemical_model")
        for name, resn in [
            (g, "LPS"),
            ("rat_Annexin_headgroup_experimental", "GSE"),
            ("PS_ideal_chemical_model", "P5S"),
        ]:
            c.color("magenta", f"{name} and resn {resn} and elem C")
            c.color("red", f"{name} and resn {resn} and elem O")
            c.color("blue", f"{name} and resn {resn} and elem N")
            c.color("orange", f"{name} and resn {resn} and elem P")
        c.color("limegreen", "elem Ca")
        c.set("sphere_scale", 0.45, "elem Ca")
        c.color("lightblue", "AXL_GAS6_experimental and chain A+B")
        c.color("orange", "AXL_GAS6_experimental and chain C+D")

        def layout(name, x, top=None, bottom=None):
            xyz = c.get_coords(name)
            centered = xyz - xyz.mean(axis=0)
            if top and bottom:
                y = c.get_coords(top).mean(axis=0) - c.get_coords(bottom).mean(axis=0)
            else:
                _, _, axes = np.linalg.svd(centered, full_matrices=False)
                y = axes[0]
            y /= np.linalg.norm(y)
            ref = np.array([0.0, 0.0, 1.0])
            if abs(ref @ y) > 0.9:
                ref = np.array([1.0, 0.0, 0.0])
            z = np.cross(ref, y)
            z /= np.linalg.norm(z)
            xx = np.cross(y, z)
            c.load_coords(centered @ np.array([xx, y, z]).T + [x, 0, 0], name)

        layout(
            "MERTK_experimental",
            -130,
            "MERTK_experimental and resi 93-186",
            "MERTK_experimental and resi 386-484",
        )
        layout(
            "GAS6_predicted",
            -55,
            "GAS6_predicted and resi 53-94",
            "GAS6_predicted and resi 477-670",
        )
        layout(
            "PROS1_predicted",
            35,
            "PROS1_predicted and resi 42-87",
            "PROS1_predicted and resi 484-666",
        )
        layout("human_Annexin_experimental", 115)
        layout("PS_ideal_chemical_model", 175)
        c.create("PROS1_LG_module_predicted", "PROS1_predicted and resi 284-676")

        scene_info = []

        def scene(name, enabled, message, focus=None, render=True):
            c.disable("all")
            for item in enabled:
                c.enable(item)
            c.deselect()
            c.orient(focus or "enabled")
            c.zoom(focus or "enabled", 5, complete=1)
            c.clip("slab", 300)
            c.scene(name, "store", message=message)
            scene_info.append({"scene": name, "message": message})
            if render:
                c.png(str(out / f"{name}.png"), width=1600, height=1000, dpi=180, ray=1)

        labels = "component_labels"
        for name, text in [
            ("MERTK_experimental", "MERTK | experimental"),
            ("GAS6_predicted", "GAS6 | predicted"),
            ("PROS1_predicted", "PROS1 | predicted"),
            ("human_Annexin_experimental", "human Annexin V"),
            ("PS_ideal_chemical_model", "PS | ideal model"),
        ]:
            xyz = c.get_coords(name)
            c.pseudoatom(
                labels, pos=[float(xyz[:, 0].mean()), float(xyz[:, 1].max() + 10), 0], label=text
            )
        c.hide("nonbonded", labels)
        scene(
            "01_separate_components",
            [
                labels,
                "MERTK_experimental",
                "GAS6_predicted",
                "PROS1_predicted",
                "human_Annexin_experimental",
                "PS_ideal_chemical_model",
            ],
            "Separate components at a common coordinate scale. Positions are for display only; this is not docking or a solved ternary complex.",
        )
        scene(
            "02_GAS6_domains",
            ["GAS6_predicted"],
            "Predicted GAS6: Gla green, EGF blue, LG1 orange, LG2 yellow-orange. No gamma-carboxylation, bound calcium or PS in this model.",
        )
        scene(
            "03_PROS1_domains",
            ["PROS1_predicted"],
            "Predicted PROS1: Gla green binds membrane; C-terminal LG module engages TAM receptors. No gamma-carboxylation, calcium or membrane represented.",
        )
        scene(
            "04_Gla_lysoPS",
            [g],
            "Experimental 1NL2: bovine prothrombin Gla, not GAS6/PROS1. Deposited CGU residues, seven calcium ions and lysoPS analogue. Other protein domains hidden.",
            f"{g} and ((polymer.protein and resi 1-46) or resn LPS or elem Ca)",
        )
        c.distance("Gla_Ca_lipid_contacts", f"{g} and elem Ca", f"{g} and resn LPS and elem O", 3.2)
        c.hide("labels", "Gla_Ca_lipid_contacts")
        scene(
            "05_Gla_contact_detail",
            [g, "Gla_Ca_lipid_contacts"],
            "Dashed Ca-O separations <=3.2 A in deposited 1NL2. Geometric contacts, not a GAS6-specific coordination model.",
            f"{g} and (resn LPS or (resn CGU within 5 of resn LPS) or (elem Ca within 5 of resn LPS))",
        )
        ann = "rat_Annexin_headgroup_experimental"
        scene(
            "06_Annexin_headgroup",
            [ann],
            "Experimental rat Annexin V 1A8A with calcium and glycerophosphoserine. This is a headgroup analogue, not a complete PS molecule or bilayer.",
        )
        c.distance(
            "Annexin_Ca_headgroup_contacts",
            f"{ann} and elem Ca",
            f"{ann} and resn GSE and elem O",
            3.2,
        )
        c.hide("labels", "Annexin_Ca_headgroup_contacts")
        c.show("sticks", f"byres ({ann} and polymer.protein within 4 of ({ann} and resn GSE))")
        scene(
            "07_Annexin_contact_detail",
            [ann, "Annexin_Ca_headgroup_contacts"],
            "Ca-O separations <=3.2 A in deposited 1A8A. Rat Annexin V; GSE lacks the PS acyl tails.",
            f"{ann} and (resn GSE or (elem Ca within 5 of resn GSE) or (polymer.protein within 5 of resn GSE))",
        )
        scene(
            "08_human_Annexin",
            ["human_Annexin_experimental"],
            "Human Annexin V/A5 1ANX, chain A, calcium-bound. No Gla domain, PS ligand or MERTK-binding module.",
        )
        scene(
            "09_receptor_binding_analogy",
            ["AXL_GAS6_experimental"],
            "Experimental 2C5D GAS6 LG module bound to AXL Ig domains. Related-receptor analogy; not a MERTK pose.",
        )
        c.show("sticks", "MERTK_experimental and polymer.protein and resi 365+367+369+460+471+472")
        c.color(
            "magenta", "MERTK_experimental and polymer.protein and resi 365+367+369+460+471+472"
        )
        scene(
            "10_MERTK_proposed_patches",
            ["MERTK_experimental"],
            "MERTK 37NT: Ig cyan/blue; FNIII teal/salmon. Proposed FNIII anchors magenta. Stalk, membrane and kinase absent; no demonstrated agonism or protection.",
        )
        scene(
            "11_PS_chemical_model",
            ["PS_ideal_chemical_model"],
            "CCD P5S ideal distearoyl PS 18:0/18:0. Chemical example, not a neutrophil-specific species, experimental bound lipid or membrane model.",
        )
        scene(
            "12_fusion_modules_separate",
            ["human_Annexin_experimental", "PROS1_LG_module_predicted"],
            "Separate candidate modules: experimental human Annexin A5 and predicted PROS1 LG region. No linker, fused coordinates, bound pose or retained activity is claimed.",
        )

        rows = []
        for name, ligand in [(g, "LPS"), (ann, "GSE")]:
            for a in c.get_model(f"{name} and elem Ca").atom:
                for b in c.get_model(f"{name} and resn {ligand} and elem O").atom:
                    d = float(np.linalg.norm(np.array(a.coord) - b.coord))
                    if d <= 3.2:
                        rows.append(
                            {
                                "object": name,
                                "Ca_residue": a.resi,
                                "lipid_residue": b.resi,
                                "lipid_atom": b.name,
                                "distance_A": round(d, 3),
                            }
                        )
        with (out / "deposited_Ca_lipid_contacts.csv").open("w") as f:
            writer = csv.DictWriter(
                f, fieldnames=["object", "Ca_residue", "lipid_residue", "lipid_atom", "distance_A"]
            )
            writer.writeheader()
            writer.writerows(rows)
        assert any(r["object"] == g for r in rows)
        assert any(r["object"] == ann for r in rows)
        c.scene("01_separate_components", "recall", animate=0)
        c.save(str(out / "PS_ligands_Annexin_MERTK.pse"))
        (out / "scenes.json").write_text(json.dumps(scene_info, indent=2) + "\n")
        print(
            f"Saved {len(scene_info)} scenes and {len(rows)} deposited Ca-lipid contacts in {out}"
        )


if __name__ == "__main__":
    main()
