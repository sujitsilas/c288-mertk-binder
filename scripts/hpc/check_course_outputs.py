#!/usr/bin/env python3
"""Audit native output geometry; never restore a distorted target to manufacture a pass."""

import argparse
import json
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser

from mertk.provenance import file_hash
from mertk.structural import rigid_fit, rmsd
from mertk.structures import read_design_chain


def structure(path):
    parser = MMCIFParser(QUIET=True) if path.suffix == ".cif" else PDBParser(QUIET=True)
    return parser.get_structure("native", path)[0]


def analyze(path, reference_path, *, context_path=None, hotspots=(365, 367, 369)):
    reference, mapping = read_design_chain(reference_path, "A")
    model = structure(path)
    chains = {}
    for chain in model:
        if any(r.id[0] == " " for r in chain):
            chains[chain.id] = read_design_chain(path, chain.id)[0]
    targets = [key for key, value in chains.items() if value.sequence == reference.sequence]
    if len(chains) != 2 or len(targets) != 1:
        raise ValueError("Need two complete protein chains and one exact reference sequence match")
    target_id = targets[0]
    binder_id = next(key for key in chains if key != target_id)
    target, binder = chains[target_id], chains[binder_id]
    rotation, translation = rigid_fit(target.backbone[:, 1], reference.backbone[:, 1])
    target_rmsd = rmsd(target.backbone[:, 1] @ rotation + translation, reference.backbone[:, 1])
    reasons = []
    if target_rmsd > 2:
        reasons.append("Target CA RMSD exceeds 2 A")
    steps = {}
    for name, chain in [("target", target), ("binder", binder)]:
        d = np.linalg.norm(np.diff(chain.backbone[:, 1], axis=0), axis=1)
        if len(d) == 0:
            raise ValueError("Need at least two residues per chain")
        steps[name] = {
            "min_A": float(d.min()),
            "max_A": float(d.max()),
            "outside_2p7_to_4p5_A": int(((d < 2.7) | (d > 4.5)).sum()),
        }
        if steps[name]["outside_2p7_to_4p5_A"]:
            reasons.append(f"{name} has abnormal adjacent CA distances")
    pairs = int(
        (
            np.linalg.norm(target.backbone[:, 1, None] - binder.backbone[None, :, 1], axis=-1) <= 8
        ).sum()
    )
    if pairs == 0:
        reasons.append("No inter-chain CA contact within 8 A")
    row = {
        "path": str(path),
        "sha256": file_hash(path),
        "target_chain": target_id,
        "binder_chain": binder_id,
        "binder_sequence": binder.sequence,
        "target_ca_rmsd_A": target_rmsd,
        "adjacent_ca": steps,
        "ca_contact_pairs": pairs,
    }
    if context_path is not None:
        context = structure(context_path)
        binder_heavy = np.array([a.coord for r in model[binder_id] for a in r if a.element != "H"])
        aligned = binder_heavy @ rotation + translation
        native = np.array(
            [a.coord for r in context["A"] if r.id[0] == " " for a in r if a.element != "H"]
        )
        glycans = np.array(
            [
                a.coord
                for a in context.get_atoms()
                if a.parent.resname in {"NAG", "MAN", "BMA"} and a.element != "H"
            ]
        )
        clash = float(np.linalg.norm(aligned[:, None] - native[None, :], axis=-1).min())
        row["intact_protein_min_heavy_atom_distance_A"] = clash
        if clash < 2.0:
            reasons.append("Clash with restored intact experimental protein context (<2 A)")
        if len(glycans):
            clearance = float(np.linalg.norm(aligned[:, None] - glycans[None, :], axis=-1).min())
            row["observed_glycan_clearance_A"] = clearance
            if clearance < 2.5:
                reasons.append("Observed glycan overlap (<2.5 A)")
        anchors = np.array([a.coord for n in hotspots for a in context["A"][n] if a.element != "H"])
        row["anchor_min_heavy_atom_distance_A"] = float(
            np.linalg.norm(aligned[:, None] - anchors[None, :], axis=-1).min()
        )
        if row["anchor_min_heavy_atom_distance_A"] > 5:
            reasons.append("No proposed anchor heavy-atom contact within 5 A")
        row["unresolved_glycan_review"] = (
            "Required; passing observed-glycan clearance cannot establish mammalian accessibility"
        )
    row.update(
        preliminary_geometry_pass=not reasons,
        rejection_reasons=reasons,
        interpretation="Raw output screening only; no affinity, AF3 nomination or functional claim",
    )
    return row


def collect(paths, reference, context=None, hotspots=(365, 367, 369)):
    rows = []
    seen = set()
    for path in paths:
        digest = file_hash(path)
        if digest in seen:
            continue
        seen.add(digest)
        try:
            rows.append(analyze(path, reference, context_path=context, hotspots=hotspots))
        except (ValueError, KeyError, IndexError) as error:
            rows.append(
                {
                    "path": str(path),
                    "sha256": digest,
                    "preliminary_geometry_pass": False,
                    "rejection_reasons": [str(error)],
                }
            )
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign", type=Path)
    parser.add_argument("--control", action="store_true")
    args = parser.parse_args()
    campaign = args.campaign.resolve()
    if args.control:
        root = campaign / "rfd4_control"
        rows = collect(
            sorted(root.glob("inference_samples/**/*_generated.cif")),
            campaign / "inputs/pdl1_control.cif",
        )
        passed = len(rows) == 1 and rows[0]["preliminary_geometry_pass"]
        (root / "quality.json").write_text(
            json.dumps({"passed": passed, "records": rows}, indent=2) + "\n"
        )
        print(
            "RFD4_CONTROL_QUALITY_PASS" if passed else "RFD4_CONTROL_QUALITY_REJECTED", flush=True
        )
        if not passed:
            raise SystemExit(42)
    else:
        manifest = json.loads((campaign / "manifest.json").read_text())
        report = {
            "status": "Preliminary geometry screen; AF3 remains pending candidate review",
            "arms": {},
        }
        for arm in ["bindcraft", "rfd4_proteina"]:
            if arm == "bindcraft":
                paths = [
                    p
                    for p in (campaign / "bindcraft/designs").rglob("*.pdb")
                    if "Trajectory" in p.parts or "Accepted" in p.parts
                ]
            else:
                paths = list((campaign / arm).glob("inference_samples/**/*_generated.cif"))
            rows = collect(
                sorted(paths),
                campaign / "inputs/target.pdb",
                campaign / "inputs/37NT.cif",
                hotspots=manifest["hotspots"],
            )
            report["arms"][arm] = {
                "records": rows,
                "files": len(rows),
                "preliminary_passes": sum(r["preliminary_geometry_pass"] for r in rows),
            }
        (campaign / "native_quality.json").write_text(json.dumps(report, indent=2) + "\n")
        print(
            json.dumps(
                {
                    arm: {k: v for k, v in data.items() if k != "records"}
                    for arm, data in report["arms"].items()
                }
            ),
            flush=True,
        )


if __name__ == "__main__":
    main()
