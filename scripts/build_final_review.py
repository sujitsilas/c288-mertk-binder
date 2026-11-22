#!/usr/bin/env python3
"""Build the offline C288 MerTK project review from the completed design runs."""

import argparse
import base64
import csv
import gzip
import html
import json
import re
import statistics
from pathlib import Path
from tempfile import NamedTemporaryFile

from Bio.PDB import MMCIFParser, PDBParser
from Bio.SeqUtils import seq1

from mertk.prediction_review import inspect_pose, safe_json
from mertk.provenance import file_hash

ROOT = Path(__file__).resolve().parents[1]
HPC = ROOT / "outputs/hpc"
FOLLOWUP = HPC / "round_20261120_mertk_pose_followup_01"
REPAIR = HPC / "round_20261121_ig1_boltz_repack_01"
IG_GUARD = HPC / "round_20261114_ig_guard_01"
RELAX = HPC / "round_20261120_ig1_relax_01"
STALK_CONFIDENCE = HPC / "round_20261117_protector_01/references/AF_Q12866.pdb"
STALK_OVERLAY_MAX_RMSD = 3.0
ASSETS = ROOT / "assets"
TEMPLATE = ASSETS / "mertk_final_review.html"
VENDOR = ASSETS / "vendor/3dmol-2.5.5"
EXPRESSION_FIGURE = ASSETS / "figures/burn_sham_mertk_axis.png"
STALK_SENSITIVITY = ASSETS / "figures/stalk_alignment_sensitivity.json"
NAMES = {
    "01_ig1_i32g_w83f": "Ig1 I32G / W83F",
    "02_ig1_q124_l129": "Ig1 Q124 · 129 aa",
    "03_ig1_q124_l114": "Ig1 Q124 · 114 aa",
    "04_gas6_positive": "GAS6",
    "05_stalk_base57_base": "Protector L57 · base",
    "06_stalk_base57_plus25": "Protector L57 · +25°",
    "07_stalk_base67_plus25": "Protector L67 · +25°",
    "08_stalk_psi25_84_plus25": "Protector L84 · +25°",
    "09_mcherry_negative": "mCherry",
}
REPAIR_NAMES = {
    "original": ("Ig1 I32G / W83F", "crystal_binder_only_202"),
    "e80q": ("Ig1 I32G / W83F + E80Q", "crystal_e80q_202"),
    "f83v": ("Ig1 I32G / W83V", "crystal_f83v_202"),
}


def read_json(path):
    return json.loads(path.read_text())


def viewer_pose(pose, alignment):
    overlay = pose["alignments"][alignment]
    return {
        "rawTarget": pose["raw_target_pdb"],
        "rawPartner": pose["raw_partner_pdb"],
        "alignedPartner": overlay["partner_pdb"],
    }


def chain_confidence(path, pose):
    parser = MMCIFParser(QUIET=True) if path.suffix == ".cif" else PDBParser(QUIET=True)
    model = parser.get_structure("prediction", path)[0]
    confidence = {}
    for label, chain in (("partnerPlddt", "partner_chain"), ("receptorPlddt", "target_chain")):
        values = [r["CA"].bfactor for r in model[pose[chain]] if r.id[0] == " "]
        if not values or not all(0.0 <= value <= 100.0 for value in values):
            raise ValueError(f"Missing or invalid Cα pLDDT values in {path}")
        confidence[label] = statistics.mean(values)
    return confidence


def stalk_confidence(canonical, sources):
    chain = PDBParser(QUIET=True).get_structure("AFDB", STALK_CONFIDENCE)[0]["A"]
    confidence = []
    for number in range(440, 531):
        residue = chain[number]
        letter = seq1(residue.resname)
        plddt = residue["CA"].bfactor
        if letter != canonical[number - 1] or not 0 < plddt <= 100:
            raise ValueError(f"AFDB stalk confidence failed sequence or value check at {number}")
        confidence.append({"number": number, "residue": letter, "plddt": plddt})
    sources[str(STALK_CONFIDENCE.relative_to(ROOT))] = file_hash(STALK_CONFIDENCE)
    return confidence


def collect_repair(cases, sources, context, canonical, receptor_sequence):
    validation = REPAIR / "validation"
    manifest_path = validation / "manifest.json"
    assessment_path = validation / "assessment.json"
    redock_path = REPAIR / "redock.json"
    jobs_path = REPAIR / "jobs.json"
    template_path = validation / "references/37NT.cif"
    manifest = read_json(manifest_path)
    template_hash = file_hash(template_path)
    if (
        template_hash != manifest["template_sha256"]
        or template_hash != sources[str((FOLLOWUP / "references/37NT.cif").relative_to(ROOT))]
    ):
        raise ValueError("Repack refolds use a different 37NT reference")
    inputs = {row["id"]: row for row in manifest["cases"]}
    redock = {row["id"]: row for row in read_json(redock_path)["cases"]}
    tested = read_json(assessment_path)["cases"]
    for path in (manifest_path, assessment_path, redock_path, jobs_path):
        sources[str(path.relative_to(ROOT))] = file_hash(path)
    sources[str(template_path.relative_to(ROOT))] = template_hash
    summary = {}
    for variant, (name, redock_id) in REPAIR_NAMES.items():
        row = redock[redock_id]
        coords = REPAIR / row["file"]
        if file_hash(coords) != row["sha256"]:
            raise ValueError(f"Restrained repair coordinates changed: {coords}")
        sources[str(coords.relative_to(ROOT))] = row["sha256"]
        summary[variant] = {
            "name": name,
            "proteinMin": row["fixed_37NT_min_A"],
            "glycanMin": row["observed_glycan_min_A"],
            "moveRmsd": row["binder_ca_rmsd_A"],
            "retained": len(row["retained_Ig1_contacts"]),
            "startingContacts": row["starting_contact_count"],
            "refolds": [],
        }
    for row in tested:
        if row["status"] == "missing":
            continue
        key = row["id"]
        variant, _ = key.rsplit("_s", 1)
        spec = inputs[key]
        directory = {"Hoffman2 CUDA": "results", "local Apple MPS": "local_results"}[row["backend"]]
        coords = validation / directory / f"boltz_results_{key}/predictions/{key}/{key}_model_0.cif"
        input_path = validation / "inputs" / f"{key}.yaml"
        if row["status"] != "complete" or file_hash(input_path) != spec["input_sha256"]:
            raise ValueError(f"Changed input or incomplete result for {key}")
        pose = inspect_pose(coords, context["Ig1"], canonical, receptor_sequence, spec["sequence"])
        if pose["source_sha256"] != row["sha256"]:
            raise ValueError(f"Refold coordinates differ from the scored result for {key}")
        all_contacts = pose["target_contacts"]
        fit = row["overlays"]["crystal"]
        name = summary[variant]["name"]
        model = {
            "id": f"repair_{key}",
            "name": f"{name} · seed {row['seed']}",
            "role": "Independent refold after restrained clash repair",
            "domain": "Ig1",
            "template": "37NT",
            "engine": "Boltz-2",
            "backend": row["backend"],
            "seed": row["seed"],
            "iptm": row["iptm"],
            "pae": row["contact_pae_A"],
            "igContacts": sum(93 <= n <= 185 for n in all_contacts),
            "stalkContacts": sum(485 <= n <= 498 for n in all_contacts),
            "cutDistances": pose["cleavage_atom_min_distances_A"],
            "proteinMin": fit["protein_min_A"],
            "glycanMin": fit["glycan_min_A"],
            "fitRmsd": fit["Ig1_fit_A"],
            "retained": len(fit["retained_contacts"]),
            "startingContacts": fit["baseline_contacts"],
            "poseRmsd": fit["partner_pose_rmsd_A"],
            "geometryPass": fit["geometry_pass"],
            "source": str(coords.relative_to(ROOT)),
            **chain_confidence(coords, pose),
            **viewer_pose(pose, "Ig1"),
        }
        cases.append(model)
        summary[variant]["refolds"].append(model["id"])
        sources[str(coords.relative_to(ROOT))] = pose["source_sha256"]
        sources[str(input_path.relative_to(ROOT))] = spec["input_sha256"]
    return {
        "variants": list(summary.values()),
        "planned": len(inputs),
        "complete": sum(len(r["refolds"]) for r in summary.values()),
        "passed": sum(r["geometryPass"] for r in cases if "geometryPass" in r),
        "job": read_json(jobs_path)["boltz_validation"],
    }


def collect_generation(sources):
    ig1_path = IG_GUARD / "campaigns/Ig1_E152/metrics.json"
    ig2_path = IG_GUARD / "campaigns/Ig2_L270/metrics.json"
    status_path = IG_GUARD / "status.json"
    relax_path = RELAX / "assessment.json"
    for path in (ig1_path, ig2_path, status_path, relax_path):
        sources[str(path.relative_to(ROOT))] = file_hash(path)
    ig1 = read_json(ig1_path)
    ig2 = read_json(ig2_path)
    relax = read_json(relax_path)["mpnn8_parent"]
    hours, minutes, seconds = map(
        int, read_json(status_path)["jobs"]["Ig1_E152_generate"]["elapsed"].split(":")
    )
    return {
        "started": ig1["attempts_started"],
        "filtered": ig1["attempts_completed"] - ig1["trajectories_scored"],
        "interrupted": ig1["attempts_censored"],
        "savedSequences": ig1["saved_sequence_proposals"],
        "acceptedSequences": ig1["accepted_variants"],
        "acceptedBackbones": ig1["accepted_trajectories"],
        "nativePasses": ig1["native_context"]["bindcraft"]["preliminary_passes"],
        "ig2Accepted": ig2["accepted_variants"],
        "ig2Sequences": ig2["saved_sequence_proposals"],
        "runtimeMinutes": hours * 60 + minutes + seconds / 60,
        "fastRelaxBefore": relax["before"]["intact_protein_min_heavy_atom_distance_A"],
        "fastRelaxAfter": relax["after"]["intact_protein_min_heavy_atom_distance_A"],
    }


def collect():
    manifest_path = FOLLOWUP / "manifest.json"
    assessment_path = FOLLOWUP / "assessment.json"
    manifest = read_json(manifest_path)
    assessment = read_json(assessment_path)
    references = FOLLOWUP / "references"
    canonical_path = references / "canonical.json"
    if file_hash(canonical_path) != manifest["canonical_source_sha256"]:
        raise ValueError("Canonical receptor sequence differs from the submitted manifest")
    canonical = read_json(canonical_path)["sequence"]["value"]
    parser = MMCIFParser(QUIET=True)
    reference_paths = {"Ig1": references / "37NT.cif", "Stalk": references / "modeled_context.cif"}
    context = {key: parser.get_structure(key, path)[0] for key, path in reference_paths.items()}
    sources = {
        str(path.relative_to(ROOT)): file_hash(path)
        for path in (manifest_path, assessment_path, canonical_path, *reference_paths.values())
    }
    reference_config = ROOT / "configs/target_reference.json"
    if (
        sources[str(reference_paths["Ig1"].relative_to(ROOT))]
        != read_json(reference_config)["source_sha256"]
    ):
        raise ValueError("Review crystal coordinates do not match the pinned 37NT source")
    sources[str(reference_config.relative_to(ROOT))] = file_hash(reference_config)
    scored = {row["id"]: row for row in assessment["cases"]}
    if set(scored) != set(NAMES) or assessment["protocol"] != manifest["protocol"]:
        raise ValueError("Unexpected final cohort or protocol")

    cases = []
    for spec in manifest["cases"]:
        key = spec["id"]
        row = scored[key]
        input_path = FOLLOWUP / "inputs" / f"{key}.yaml"
        if file_hash(input_path) != spec["input_sha256"] or row["status"] != "complete":
            raise ValueError(f"Changed input or incomplete result for {key}")
        coords = FOLLOWUP / row["structure"]
        domain = "Stalk" if key[3:].startswith("stalk_") else "Ig1"
        pose = inspect_pose(
            coords, context[domain], canonical, manifest["receptor_sequence"], spec["sequence"]
        )
        if pose["source_sha256"] != row["source_sha256"]:
            raise ValueError(f"Coordinates differ from the scored result for {key}")
        audit = row["overlays"]["modeled_stalk" if domain == "Stalk" else "experimental_Ig1"]
        stalk_contacts = len(row["stalk_contact_residues_4p5A"])
        cases.append(
            {
                "id": key,
                "name": NAMES[key],
                "role": row["role"],
                "domain": domain,
                "template": spec["target_template"],
                "engine": "Boltz-2",
                "iptm": row["iptm"],
                "pairIptm": [row["pair_chains_iptm"][a][b] for a, b in (("0", "1"), ("1", "0"))],
                "pae": row["mean_bidirectional_contact_PAE_A"],
                "igContacts": len(row["Ig1_contact_residues_4p5A"]),
                "stalkContacts": stalk_contacts,
                "cutDistances": row["P490_C_and_S491_N_min_A"],
                "proteinMin": audit["protein_min_A"],
                "glycanMin": audit["glycan_min_A"],
                "fitRmsd": audit["target_ca_rmsd_A"],
                "overlayReliable": domain != "Stalk"
                or (audit["target_ca_rmsd_A"] <= STALK_OVERLAY_MAX_RMSD and stalk_contacts > 0),
                "retained": (
                    len(row["retained_source_pose_contacts"])
                    if "retained_source_pose_contacts" in row
                    else None
                ),
                "startingContacts": row.get("source_pose_contact_count"),
                "source": str(coords.relative_to(ROOT)),
                **chain_confidence(coords, pose),
                **viewer_pose(pose, domain),
            }
        )
        sources[str(coords.relative_to(ROOT))] = pose["source_sha256"]
        sources[str(input_path.relative_to(ROOT))] = spec["input_sha256"]

    ig1 = IG_GUARD / "campaigns/Ig1_E152/bindcraft/designs"
    statistics_path = ig1 / "final_design_stats.csv"
    with statistics_path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            variant = row["Design"].split("_")[-1]
            if variant not in {"mpnn8", "mpnn12"}:
                continue
            (coords,) = (ig1 / "Accepted").glob(f"{row['Design']}_model*.pdb")
            pose = inspect_pose(coords, context["Ig1"], canonical)
            audit = pose["alignments"]["Ig1"]
            cases.append(
                {
                    "id": variant,
                    "name": f"Ig1 {variant} · original",
                    "role": "Original E152 BindCraft design",
                    "domain": "Ig1",
                    "engine": "BindCraft / AF2",
                    "iptm": float(row["Average_i_pTM"]),
                    "proteinMin": audit["protein_min_A"],
                    "glycanMin": audit["glycan_min_A"],
                    "fitRmsd": audit["target_ca_rmsd_A"],
                    "source": str(coords.relative_to(ROOT)),
                    **chain_confidence(coords, pose),
                    **viewer_pose(pose, "Ig1"),
                }
            )
            sources[str(coords.relative_to(ROOT))] = pose["source_sha256"]

    local = HPC / "round_20261120_ig1_local_geometry_02"
    geometry_path = local / "geometry.json"
    refold_path = local / "independent.json"
    repair_dir = HPC / "round_20261119_proposal_01/results/repairs"
    repair_paths = sorted(repair_dir.glob("*/*/metrics.json"))
    if not repair_paths:
        raise ValueError("Missing independent RFD4 repair scores")
    stalk_path = FOLLOWUP / "stalk_geometry.json"
    job_path = FOLLOWUP / "jobs.json"
    for path in (statistics_path, geometry_path, refold_path, stalk_path, job_path, *repair_paths):
        sources[str(path.relative_to(ROOT))] = file_hash(path)
    sources[str(STALK_SENSITIVITY.relative_to(ROOT))] = file_hash(STALK_SENSITIVITY)
    staged = read_json(geometry_path)
    refolds = read_json(refold_path)
    stalk = read_json(stalk_path)
    workflow = {
        "localPassed": sum(row["pose_pass"] for row in staged.values()),
        "localTotal": len(staged),
        "refoldPassed": sum(row["pose_pass"] for row in refolds.values()),
        "refoldTotal": len(refolds),
        "rfdPassed": sum(
            read_json(p)["native_context"]["preliminary_geometry_pass"] for p in repair_paths
        ),
        "rfdTotal": len(repair_paths),
        "stalkTotal": len(stalk),
        "stalkPassed": sum(all(c["pass"] for c in row["conformers"]) for row in stalk),
        "boltzTotal": len(scored),
    }
    return {
        "cases": cases,
        "generation": collect_generation(sources),
        "repair": collect_repair(cases, sources, context, canonical, manifest["receptor_sequence"]),
        "references": {key: path.read_text() for key, path in reference_paths.items()},
        "stalkConfidence": stalk_confidence(canonical, sources),
        "stalkSensitivity": read_json(STALK_SENSITIVITY),
        "stalkOverlayMaxRmsd": STALK_OVERLAY_MAX_RMSD,
        "workflow": workflow,
        "protocol": manifest["protocol"],
        "range": manifest["canonical_range"],
        "jobs": read_json(job_path),
    }, sources


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--public", action="store_true", help="Build a single page for GitHub Pages"
    )
    args = parser.parse_args()
    if args.public and args.output is None:
        parser.error("--public requires --output")
    output = args.output or ROOT / "outputs/review/mertk_protector_predictions.html"
    data, sources = collect()
    library = VENDOR / "3Dmol-min.js"
    digest = file_hash(library)
    license_path = VENDOR / "LICENSE"
    vendor_files = read_json(VENDOR / "provenance.json")["files"]
    if (
        digest != vendor_files[library.name]["sha256"]
        or file_hash(license_path) != vendor_files[license_path.name]["sha256"]
    ):
        raise ValueError("Pinned 3Dmol.js library or license hash changed")
    figure = EXPRESSION_FIGURE.read_bytes()
    sources[str(EXPRESSION_FIGURE.relative_to(ROOT))] = file_hash(EXPRESSION_FIGURE)
    template = TEMPLATE.read_text()
    if args.public:
        template = re.sub(r"<!-- local-only -->.*?<!-- /local-only -->", "", template, flags=re.S)
        for case in data["cases"]:
            del case["source"]
    compressed_data = gzip.compress(safe_json(data).encode("utf-8"), mtime=0)
    page = (
        template.replace("__DATA__", base64.b64encode(compressed_data).decode("ascii"))
        .replace("__LIBRARY__", library.read_text())
        .replace("__LICENSE__", html.escape(license_path.read_text()))
        .replace("__EXPRESSION_FIGURE__", base64.b64encode(figure).decode("ascii"))
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile(mode="w", encoding="utf-8", dir=output.parent, delete=False) as stream:
        stream.write(page)
        staged = Path(stream.name)
    staged.chmod(0o644)
    staged.replace(output)
    if args.public:
        print(f"{output} · {len(data['cases'])} structures")
        return
    summary = {
        "workflow": data["workflow"],
        "generation": data["generation"],
        "repair": data["repair"],
        "stalkConfidence": data["stalkConfidence"],
        "stalkSensitivity": data["stalkSensitivity"],
        "stalkOverlayMaxRmsd": data["stalkOverlayMaxRmsd"],
        "jobs": data["jobs"],
        "cases": [
            {
                k: v
                for k, v in case.items()
                if k not in {"rawTarget", "rawPartner", "alignedPartner"}
            }
            for case in data["cases"]
        ],
        "sources": sources,
        "template_sha256": file_hash(TEMPLATE),
        "library_sha256": digest,
        "html_sha256": file_hash(output),
    }
    output.with_suffix(".json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"{output} · {len(data['cases'])} structures")


if __name__ == "__main__":
    main()
