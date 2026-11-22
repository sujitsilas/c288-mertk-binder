"""Summarize completed Boltz protein complexes, retaining missing or incomplete cases."""

import argparse
import json
import sys
from inspect import getsourcefile
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser

from mertk.prediction_review import inspect_pose
from mertk.provenance import file_hash, source_revision


def score(suite, case, target, canonical, native, modeled):
    name = case["id"]
    directory = suite / "results" / f"boltz_results_{name}" / "predictions" / name
    coords = directory / f"{name}_model_0.cif"
    confidence = directory / f"confidence_{name}_model_0.json"
    pae_path = directory / f"pae_{name}_model_0.npz"
    missing = [p.name for p in (coords, confidence, pae_path) if not p.is_file()]
    if missing:
        return {"status": "incomplete", "missing": missing}
    pose = inspect_pose(coords, native, canonical, target, case["sequence"])
    stalk_pose = inspect_pose(coords, modeled, canonical, target, case["sequence"])
    model = MMCIFParser(QUIET=True).get_structure("boltz", coords)[0]
    if [c.id for c in model] != ["A", "B"] or pose["target_chain"] != "A":
        raise ValueError("Expected receptor A followed by partner B for the PAE token order")
    receptor = list(model["A"])
    partner = list(model["B"])
    x = np.array([r["CA"].coord for r in receptor])
    y = np.array([r["CA"].coord for r in partner])
    i, j = np.where(np.linalg.norm(x[:, None] - y, axis=-1) <= 10)
    with np.load(pae_path) as errors:
        pae = errors["pae"]
        if pae.shape != (len(x) + len(y),) * 2:
            raise ValueError("PAE dimensions differ from the two protein chains")
        contact_pae = (
            float((pae[i, j + len(x)] + pae[j + len(x), i]).mean() / 2) if len(i) else None
        )
    ligand_xyz = np.array([a.coord for a in model["B"].get_atoms() if a.element != "H"])
    target_start = canonical.index(target) + 1
    anchor = np.array([a.coord for a in receptor[152 - target_start] if a.element != "H"])
    confidence_data = json.loads(confidence.read_text())
    contacts = pose["target_contacts"]
    overlays = {}
    for label, fit in [
        ("experimental_Ig1", pose["alignments"]["Ig1"]),
        ("modeled_stalk", stalk_pose["alignments"]["Stalk"]),
    ]:
        overlays[label] = {
            key: fit[key] for key in ["protein_min_A", "glycan_min_A", "target_ca_rmsd_A"]
        }
    return {
        "status": "complete",
        "structure": str(coords.relative_to(suite)),
        "source_sha256": file_hash(coords),
        "iptm": confidence_data["iptm"],
        "pair_chains_iptm": confidence_data["pair_chains_iptm"],
        "receptor_partner_CA_pairs_within_10A": len(i),
        "mean_bidirectional_contact_PAE_A": contact_pae,
        "Ig1_contact_residues_4p5A": [n for n in contacts if 93 <= n <= 185],
        "stalk_contact_residues_4p5A": [n for n in contacts if 485 <= n <= 498],
        "E152_min_heavy_atom_A": float(np.linalg.norm(ligand_xyz[:, None] - anchor, axis=-1).min()),
        "P490_C_and_S491_N_min_A": stalk_pose["cleavage_atom_min_distances_A"],
        "overlays": overlays,
    }


def main(suite):
    manifest = json.loads((suite / "manifest.json").read_text())
    references = suite / "references"
    canonical_file = references / "canonical.json"
    if file_hash(canonical_file) != manifest["canonical_source_sha256"]:
        raise ValueError("Unexpected canonical reference")
    canonical = json.loads(canonical_file.read_text())["sequence"]["value"]
    parser = MMCIFParser(QUIET=True)
    native = parser.get_structure("37NT", references / "37NT.cif")[0]
    modeled = parser.get_structure("stalk", references / "modeled_context.cif")[0]
    rows = []
    for case in manifest["cases"]:
        if file_hash(suite / "inputs" / (case["id"] + ".yaml")) != case["input_sha256"]:
            raise ValueError("Prediction input differs from the submitted manifest")
        try:
            row = score(suite, case, manifest["receptor_sequence"], canonical, native, modeled)
            if row["status"] == "complete" and "reference_contact_residues_4p5A" in case:
                if case["role"].startswith("Aim 1"):
                    site_key = "Ig1_contact_residues_4p5A"
                elif case["role"].startswith("Aim 2"):
                    site_key = "stalk_contact_residues_4p5A"
                else:
                    raise ValueError(f"Reference contact region is unspecified for {case['role']}")
                baseline = set(case["reference_contact_residues_4p5A"])
                row["source_pose_contact_count"] = len(baseline)
                row["retained_source_pose_contacts"] = sorted(baseline.intersection(row[site_key]))
        except (KeyError, ValueError, OSError, IndexError) as error:
            row = {"status": "unscored", "error": str(error)}
        rows.append(
            {
                "id": case["id"],
                "role": case["role"],
                "target_template": case.get("target_template"),
                **row,
            }
        )
    summary = {
        "protocol": manifest["protocol"],
        "scorer": {
            "sha256": file_hash(__file__),
            "pose_inspection_sha256": file_hash(getsourcefile(inspect_pose)),
            **source_revision(__file__),
        },
        "cases": rows,
    }
    path = suite / "assessment.json"
    path.with_suffix(".tmp").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    path.with_suffix(".tmp").replace(path)
    print([(r["id"], r["status"]) for r in rows], flush=True)
    return int(any(r["status"] != "complete" for r in rows))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("suite", type=Path)
    sys.exit(main(parser.parse_args().suite))
