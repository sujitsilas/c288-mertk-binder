#!/usr/bin/env python3
"""Prepare receptor-only Boltz replicates and score Ig1 pose recovery after repacking."""

import argparse
import copy
import json
from io import StringIO
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser
from repack_ig1_boltz import LYSINES, SEEDS, coords, distance, heavy, proteins, save_json

from mertk.prediction_review import inspect_pose
from mertk.provenance import file_hash, source_revision
from mertk.structural import rmsd


def rank(row):
    return (
        -len(row["retained_Ig1_contacts"]),
        -row["closest_protein_A"],
        -row["observed_glycan_min_A"],
    )


def prepare(suite, original):
    metadata = json.loads((suite / "inputs.json").read_text())
    rows = json.loads((suite / "repack.json").read_text())["cases"]
    passing = [row for row in rows if row["geometry_pass"]]
    unchanged = [row for row in passing if row["protocol"] == "both_sides"]
    variants = [
        row
        for row in passing
        if row["binder_sequence"] != metadata["expected_binder"] and row["fixed_37NT_min_A"] >= 2
    ]
    selected = []
    if unchanged:
        selected.append(("original", min(unchanged, key=rank)))
    # Verify one substitution for each clash, since their fixed-backbone geometries tied.
    for label in ("e80q", "f83v"):
        eligible = [row for row in variants if row["protocol"] == label]
        if eligible:
            selected.append((label, min(eligible, key=rank)))
    if not selected:
        raise ValueError("No geometry-passing repacked pose to validate")
    # Boltz receives receptor template and sequence only. No repacked coordinates constrain docking.
    original_yaml = original.read_text()
    binder_line = f"sequence: {metadata['expected_binder']}"
    if original_yaml.count(binder_line) != 1 or "cif: references/37NT.cif" not in original_yaml:
        raise ValueError("Starting Boltz YAML changed")
    work = suite / "validation"
    inp = work / "inputs"
    inp.mkdir(parents=True, exist_ok=False)
    (work / "references").mkdir()
    (work / "references/37NT.cif").write_bytes((suite / "37NT.cif").read_bytes())
    cases = []
    for label, parent in selected:
        for seed in SEEDS:
            name = f"{label}_s{seed}"
            path = inp / f"{name}.yaml"
            path.write_text(
                original_yaml.replace(binder_line, f"sequence: {parent['binder_sequence']}")
            )
            cases.append(
                {
                    "id": name,
                    "seed": seed,
                    "sequence": parent["binder_sequence"],
                    "parent": parent["file"],
                    "parent_sha256": parent["sha256"],
                    "input_sha256": file_hash(path),
                }
            )
    save_json(
        work / "manifest.json",
        {
            "cases": cases,
            "source_input_sha256": file_hash(original),
            "template_sha256": file_hash(work / "references/37NT.cif"),
            "starting_contacts": metadata["starting_contacts"],
            "protocol": {
                "model": "Boltz-2",
                "samples": 1,
                "recycles": 3,
                "sampling_steps": 200,
                "receptor_template_only": True,
                **source_revision(__file__),
                "workflow_sha256": file_hash(__file__),
            },
        },
    )


def overlay_result(fit, reference, baseline, source_binder):
    parser = PDBParser(QUIET=True)
    aligned = parser.get_structure("aligned", StringIO(fit["partner_pdb"]))[0]["B"]
    xyz = coords(heavy(aligned.get_atoms()))
    retained = [n for n in baseline if distance(coords(heavy(reference["A"][n])), xyz).min() <= 4.5]
    partner_rmsd = rmsd(
        coords(r["CA"] for r in proteins(aligned)), coords(r["CA"] for r in proteins(source_binder))
    )
    return {
        "protein_min_A": fit["protein_min_A"],
        "glycan_min_A": fit["glycan_min_A"],
        "Ig1_fit_A": fit["target_ca_rmsd_A"],
        "partner_pose_rmsd_A": partner_rmsd,
        "retained_contacts": retained,
        "baseline_contacts": len(baseline),
        "geometry_pass": fit["protein_min_A"] >= 2
        and fit["glycan_min_A"] >= 2.5
        and len(retained) >= 0.7 * len(baseline)
        and partner_rmsd <= 2.5,
    }


def score(suite):
    work = suite / "validation"
    manifest = json.loads((work / "manifest.json").read_text())
    canonical = json.loads((suite / "canonical.json").read_text())["sequence"]["value"]
    template = work / "references/37NT.cif"
    if file_hash(template) != manifest["template_sha256"]:
        raise ValueError("Frozen validation template changed")
    native = MMCIFParser(QUIET=True).get_structure("native", template)[0]
    parser = PDBParser(QUIET=True)
    result = []
    for case in manifest["cases"]:
        name = case["id"]
        for backend, directory in [
            ("Hoffman2 CUDA", "results"),
            ("local Apple MPS", "local_results"),
        ]:
            pred = work / directory / f"boltz_results_{name}" / "predictions" / name
            path = pred / f"{name}_model_0.cif"
            confidence = pred / f"confidence_{name}_model_0.json"
            errors = pred / f"pae_{name}_model_0.npz"
            if all(p.is_file() for p in (path, confidence, errors)):
                break
        else:
            result.append({"id": name, "status": "missing"})
            continue
        parent = suite / case["parent"]
        yaml = work / "inputs" / f"{name}.yaml"
        if (file_hash(parent), file_hash(yaml)) != (case["parent_sha256"], case["input_sha256"]):
            raise ValueError(f"Inputs changed for {name}")
        frozen = parser.get_structure(name, parent)[0]
        repacked = copy.deepcopy(native)
        for n in LYSINES:
            for atom in repacked["A"][n]:
                atom.coord = frozen["A"][n][atom.name].coord.copy()
        overlays = {}
        for label, reference in (("crystal", native), ("repacked_lysines", repacked)):
            pose = inspect_pose(path, reference, canonical, canonical[20:505], case["sequence"])
            overlays[label] = overlay_result(
                pose["alignments"]["Ig1"], reference, manifest["starting_contacts"], frozen["B"]
            )
        raw = MMCIFParser(QUIET=True).get_structure("prediction", path)[0]
        a, b = [coords(r["CA"] for r in proteins(raw[chain])) for chain in ("A", "B")]
        i, j = np.where(distance(a, b) <= 10)
        with np.load(errors) as arrays:
            pae = arrays["pae"]
            if pae.shape != (len(a) + len(b),) * 2:
                raise ValueError(f"Unexpected PAE shape for {name}")
            contact_pae = (
                float((pae[i, j + len(a)] + pae[j + len(a), i]).mean() / 2) if len(i) else None
            )
        result.append(
            {
                "id": name,
                "status": "complete",
                "backend": backend,
                "seed": case["seed"],
                "iptm": json.loads(confidence.read_text())["iptm"],
                "contact_pae_A": contact_pae,
                "raw_protein_min_A": float(
                    distance(
                        coords(heavy(raw["A"].get_atoms())), coords(heavy(raw["B"].get_atoms()))
                    ).min()
                ),
                "sha256": file_hash(path),
                "overlays": overlays,
            }
        )
    save_json(work / "assessment.json", {"cases": result})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["prepare", "score"])
    parser.add_argument("suite", type=Path)
    parser.add_argument("--original", type=Path)
    args = parser.parse_args()
    if args.stage == "prepare":
        prepare(args.suite, args.original)
    else:
        score(args.suite)


if __name__ == "__main__":
    main()
