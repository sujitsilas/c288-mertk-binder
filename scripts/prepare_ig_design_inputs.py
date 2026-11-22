#!/usr/bin/env python3
"""Prepare Ig-only target coordinates and explicit glycan exclusion review inputs.

This is a CPU target-preparation step, not an adapter for a design engine. A residue
mask is a heuristic; recorded carbohydrate coordinates must also constrain whole poses.
"""

import argparse
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from Bio.PDB import PDBIO, MMCIFParser, PDBParser, Select

from mertk.provenance import file_hash
from mertk.structures import read_design_chain
from mertk.target_policy import require_ig_target

SUGARS = {"NAG", "BMA", "MAN", "FUC", "GAL", "NDG", "A2G", "SIA"}


class IgPair(Select):
    def accept_chain(self, chain):
        return chain.id == "A"

    def accept_residue(self, residue):
        return residue.id[0] == " " and 93 <= residue.id[1] <= 285


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("Use a new output directory; prior inputs remain immutable")
    reference = Path("data/references/37NT")
    provenance = json.loads((reference / "provenance.json").read_text())
    for name, key in [("target.pdb", "target_sha256"), ("37NT.cif", "source_sha256")]:
        if file_hash(reference / name) != provenance[key]:
            raise ValueError("Experimental reference checksum changed")
    structure = PDBParser(QUIET=True).get_structure("native", reference / "target.pdb")
    native = structure[0]["A"]
    glyco = MMCIFParser(QUIET=True).get_structure("glyco", reference / "37NT.cif")[0]
    glycan_atoms = [
        {
            "chain": a.parent.parent.id,
            "residue": a.parent.id[1],
            "resname": a.parent.resname,
            "atom": a.name,
            "element": a.element,
            "xyz_A": a.coord.tolist(),
        }
        for r in glyco.get_residues()
        if r.resname in SUGARS
        for a in r
        if a.element != "H"
    ]
    xyz = np.array([a["xyz_A"] for a in glycan_atoms])
    full, full_map = read_design_chain(reference / "target.pdb", "A")
    sequence = {r["author_number"]: aa for r, aa in zip(full_map, full.sequence)}
    sequons = [
        n
        for n in sorted(sequence)
        if sequence[n] == "N"
        and sequence.get(n + 1) not in (None, "P")
        and sequence.get(n + 2) in ("S", "T")
    ]
    epitopes = json.loads(Path("configs/course_epitopes.json").read_text())
    scaffold = json.loads(Path("configs/ig_scaffold_design.json").read_text())
    sites = {}
    for key in scaffold["target_epitopes"]:
        site = epitopes[key]
        domain = require_ig_target(key, site["hotspots"])
        rows = []
        for n in range(93, 286):
            atom_dist = float(
                np.linalg.norm(
                    np.array([a.coord for a in native[n] if a.element != "H"])[:, None] - xyz[None],
                    axis=-1,
                ).min()
            )
            distance, nearest = min(
                (float(np.linalg.norm(native[n]["CA"].coord - native[s]["CA"].coord)), s)
                for s in sequons
            )
            rows.append(
                {
                    "residue": n,
                    "nearest_sequon": nearest,
                    "sequon_ca_distance_A": round(distance, 3),
                    "observed_glycan_heavy_distance_A": round(atom_dist, 3),
                    "avoid": distance < 15.0 or atom_dist < 6.0,
                }
            )
        coldspots = [r["residue"] for r in rows if r["avoid"]]
        overlap = sorted(set(coldspots) & set(site["hotspots"]))
        if overlap:
            raise ValueError(f"{key} hotspots conflict with conservative glycan mask: {overlap}")
        centre = np.mean([native[n]["CA"].coord for n in site["hotspots"]], axis=0)
        clearance = min(
            float(np.linalg.norm(a.coord - centre))
            for r in native
            if r.id[0] == " " and r.id[1] > 285
            for a in r
            if a.element != "H"
        )
        if clearance < 20:
            raise ValueError(f"{key} requires a larger receptor crop")
        sites[key] = {
            **site,
            "domain": domain,
            "crop_clearance_A": clearance,
            "coldspot_residues": coldspots,
            "residue_mask_audit": rows,
            "generator_residue_map": {str(n - 92): n for n in range(93, 286)},
            "generation_ready": False,
        }
    args.output.mkdir(parents=True)
    writer = PDBIO()
    writer.set_structure(structure)
    writer.save(str(args.output / "ig_pair.pdb"), select=IgPair())
    crop, mapping = read_design_chain(args.output / "ig_pair.pdb", "A")
    if [r["author_number"] for r in mapping] != list(range(93, 286)):
        raise ValueError("Crop numbering changed")
    if crop.sequence != full.sequence[:193] or not np.allclose(crop.backbone, full.backbone[:193]):
        raise ValueError("Native crop sequence or coordinates changed")
    shutil.copy2(reference / "37NT.cif", args.output / "full_native_glycosylated.cif")
    shutil.copy2(Path(__file__), args.output / "preparation_script.py")
    for name in [
        "active_target_policy.json",
        "ig_scaffold_design.json",
        "project_aims.json",
        "mertk_cleavage_protector.json",
    ]:
        shutil.copy2(Path("configs") / name, args.output / name)
    (args.output / "glycan_exclusion_review.json").write_text(
        json.dumps(
            {
                "status": "Prepared review constraints; not consumed by the current generator",
                "coordinates": "Native 37NT author coordinate frame; do not compare to renumbered or moved poses without sequence mapping and alignment",
                "mask_rule": "Target residue CA <15 A from any sequon CA OR any heavy atom <6 A from an observed glycan heavy atom",
                "limitations": "Conservative heuristic, not a complete mammalian glycan ensemble or physical energy model. Coldspots alone do not exclude carbohydrate volume.",
                "sequons": [
                    {
                        "residue": n,
                        "ca_xyz_A": native[n]["CA"].coord.tolist(),
                        "nd2_xyz_A": native[n]["ND2"].coord.tolist(),
                    }
                    for n in sequons
                ],
                "observed_glycan_atoms": glycan_atoms,
                "sites": sites,
            },
            indent=2,
        )
        + "\n"
    )
    record = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "purpose": "Ig1/Ig2-only target and scaffold preparation; no generation or submission",
        "epitopes": list(sites),
        "generation_ready": False,
        "source_sha256": {
            str(p): file_hash(p)
            for p in [
                reference / "target.pdb",
                reference / "37NT.cif",
                Path(__file__),
                Path("configs/course_epitopes.json"),
                Path("configs/ig_scaffold_design.json"),
                Path("configs/project_aims.json"),
                Path("configs/mertk_cleavage_protector.json"),
            ]
        },
        "output_sha256": {p.name: file_hash(p) for p in args.output.iterdir() if p.is_file()},
    }
    (args.output / "manifest.json").write_text(json.dumps(record, indent=2) + "\n")
    print(
        json.dumps(
            {
                "output": str(args.output),
                "epitopes": list(sites),
                "sequons": len(sequons),
                "glycan_atoms": len(glycan_atoms),
                "generation_ready": False,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
