#!/usr/bin/env python3
"""Freeze configured BindCraft experiments and their workflow; never overwrite an old run."""

import argparse
import json
import re
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath

import numpy as np
from Bio.PDB import PDBIO, PDBParser, Select

from mertk.glycan_guard import prepare_constraints
from mertk.provenance import file_hash
from mertk.structures import read_design_chain
from mertk.target_policy import require_ig_target


class Crop(Select):
    def __init__(self, start, end):
        self.start, self.end = start, end

    def accept_chain(self, chain):
        return chain.id == "A"

    def accept_residue(self, residue):
        return residue.id[0] == " " and self.start <= residue.id[1] <= self.end


def crop_clearance(chain, limits, hotspots):
    centre = np.mean([chain[n]["CA"].coord for n in hotspots], axis=0)
    excluded = sorted(
        (float(min(np.linalg.norm(a.coord - centre) for a in r if a.element != "H")), r.id[1])
        for r in chain
        if r.id[0] == " " and not limits[0] <= r.id[1] <= limits[1]
    )
    if excluded and excluded[0][0] < 20:
        raise ValueError(f"Crop excludes residues within 20 A of hotspot centroid: {excluded[:5]}")
    return {
        "nearest_excluded_heavy_atom_A": excluded[0][0] if excluded else None,
        "nearest_excluded_residue": excluded[0][1] if excluded else None,
        "minimum_required_A": 20,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", required=True)
    parser.add_argument("--remote-project", required=True)
    parser.add_argument("--config", type=Path, default=Path("configs/bindcraft_round.json"))
    parser.add_argument("--output-root", type=Path, default=Path("outputs/hpc"))
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9_-]+", args.campaign):
        parser.error("Invalid campaign identifier")
    project = PurePosixPath(args.remote_project)
    if not project.is_absolute():
        parser.error("Remote project must be absolute")
    out = args.output_root / args.campaign
    if out.exists():
        raise FileExistsError("Existing campaign is immutable; choose a fresh identifier")
    config = json.loads(args.config.read_text())
    epitopes = json.loads(Path("configs/course_epitopes.json").read_text())
    for key, spec in config["experiments"].items():
        if key == "pdl1_control" and spec["epitope"] is None:
            continue
        epitope = epitopes[spec["epitope"]]
        require_ig_target(spec["epitope"], epitope["hotspots"])
    tools = json.loads(Path("configs/external_tools.json").read_text())
    engine = Path(tools["bindcraft"]["directory"])
    revision = subprocess.check_output(
        ["git", "-C", str(engine), "rev-parse", "HEAD"], text=True
    ).strip()
    if revision != tools["bindcraft"]["revision"]:
        raise ValueError("BindCraft checkout is not pinned")
    reference = Path("data/references/37NT")
    provenance = json.loads((reference / "provenance.json").read_text())
    for name, key in [("37NT.cif", "source_sha256"), ("target.pdb", "target_sha256")]:
        if file_hash(reference / name) != provenance[key]:
            raise ValueError("Reference hash changed")
    structure = PDBParser(QUIET=True).get_structure("mertk", reference / "target.pdb")
    full, full_map = read_design_chain(reference / "target.pdb", "A")
    # Validate every crop before creating a partial campaign directory.
    clearances = {
        key: crop_clearance(
            structure[0]["A"],
            epitopes[spec["epitope"]]["crop"],
            epitopes[spec["epitope"]]["hotspots"],
        )
        for key, spec in config["experiments"].items()
        if spec["epitope"]
    }
    remote = project / "runs" / args.campaign
    workflow = out / "workflow"
    workflow.mkdir(parents=True)
    for name in [
        "run_bounded_bindcraft.py",
        "run_bindcraft_round.sh",
        "bindcraft_round.py",
        "check_course_outputs.py",
        "prepare_bindcraft_round.py",
        "submit_bindcraft_round.py",
    ]:
        shutil.copy2(Path("scripts/hpc") / name, workflow / name)
    if config.get("workflow") == "ig_followup_v1":
        for name in [
            "check_glycan_runtime.py",
            "run_ig_followup.sh",
            "ig_followup.py",
            "submit_ig_followup.py",
        ]:
            shutil.copy2(Path("scripts/hpc") / name, workflow / name)
        control = Path(config["reused_positive_control"])
        proof = json.loads(control.read_text())
        if proof.get("passed") is not True:
            raise ValueError("Reused positive control did not pass")
        shutil.copy2(control, workflow / "positive_control_evidence.json")
        comparisons = out / "comparisons"
        comparisons.mkdir()
        source = Path(config["controls_manifest"])
        controls = json.loads(source.read_text())
        for case in controls["cases"]:
            original = source.parent / case["input"]
            if file_hash(original) != case["input_sha256"]:
                raise ValueError("Comparator sequence input changed")
            destination = comparisons / case["input"]
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(original, destination)
            case["status"] = "Prepared for queued AF3 pilot"
        controls["status"] = (
            "Matched three-seed AF3 pilot; user requested controls before qualified parents"
        )
        controls["protocol"]["state"] = "locked_pilot"
        controls["protocol"]["flash_attention_implementation"] = "triton"
        controls["release_requirements"][0] = (
            "User authorized comparator inference now; any later binder comparison must use this same protocol."
        )
        (comparisons / "manifest.json").write_text(json.dumps(controls, indent=2) + "\n")
        shutil.copy2("configs/external_tools.json", workflow / "external_tools.json")
    shutil.copytree("src/mertk", workflow / "mertk", ignore=shutil.ignore_patterns("__pycache__"))
    shutil.copy2(args.config, workflow / "round_config.json")
    shutil.copy2("configs/course_epitopes.json", workflow / "course_epitopes.json")
    shutil.copy2(
        "configs/active_target_policy.json", workflow / "target_policy_at_preparation.json"
    )
    shutil.copy2("docs/efferocytosis_design.md", workflow / "proposal_at_preparation.md")
    for key, spec in config["experiments"].items():
        campaign = out / "campaigns" / key
        inputs = campaign / "inputs"
        (inputs / "bindcraft").mkdir(parents=True)
        experiment_remote = remote / "campaigns" / key
        epitope = epitopes[spec["epitope"]] if spec["epitope"] else None
        if epitope:
            writer = PDBIO()
            writer.set_structure(structure)
            writer.save(str(inputs / "target.pdb"), select=Crop(*epitope["crop"]))
            crop, mapping = read_design_chain(inputs / "target.pdb", "A")
            numbers = [r["author_number"] for r in mapping]
            expected = list(range(epitope["crop"][0], epitope["crop"][1] + 1))
            indices = [i for i, r in enumerate(full_map) if r["author_number"] in expected]
            if numbers != expected or crop.sequence != "".join(full.sequence[i] for i in indices):
                raise ValueError("Crop residue mapping or sequence changed")
            if not np.allclose(crop.backbone, full.backbone[indices]):
                raise ValueError("Crop coordinates changed")
            shutil.copy2(reference / "37NT.cif", inputs / "37NT.cif")
            shutil.copy2(reference / "provenance.json", inputs / "reference_provenance.json")
            (inputs / "residue_map.json").write_text(json.dumps(mapping, indent=2) + "\n")
            if config.get("workflow") == "ig_followup_v1":
                constraints = prepare_constraints(inputs / "target.pdb", inputs / "37NT.cif")
                (inputs / "glycan_constraints.json").write_text(
                    json.dumps(constraints, indent=2) + "\n"
                )
            settings = {
                "chains": "A",
                "target_hotspot_residues": ",".join(f"A{n}" for n in epitope["hotspots"]),
            }
        else:
            shutil.copy2(engine / "example/PDL1.pdb", inputs / "target.pdb")
            settings = json.loads((engine / "settings_target/PDL1.json").read_text())
        settings.update(
            design_path=str(experiment_remote / "bindcraft/designs"),
            binder_name=key,
            starting_pdb=str(experiment_remote / "inputs/target.pdb"),
            lengths=spec["lengths"],
            number_of_final_designs=spec["final_designs"],
        )
        advanced = json.loads((engine / "settings_advanced" / spec["advanced_preset"]).read_text())
        advanced.update(
            max_trajectories=spec["max_attempts"],
            num_seqs=spec["num_seqs"],
            max_mpnn_sequences=spec["max_mpnn_sequences"],
            save_design_animations=False,
            save_design_trajectory_plots=True,
            save_mpnn_fasta=True,
            remove_unrelaxed_trajectory=False,
            remove_unrelaxed_complex=False,
            remove_binder_monomer=False,
            af_params_dir=str(project / "external/BindCraft"),
        )
        for name, data in [("target", settings), ("advanced", advanced)]:
            (inputs / "bindcraft" / f"{name}.json").write_text(json.dumps(data, indent=2) + "\n")
        shutil.copy2(
            engine / "settings_filters/default_filters.json", inputs / "bindcraft/filters.json"
        )
        record = {
            "campaign_id": f"{args.campaign}_{key}",
            "remote_project": str(project),
            "remote_campaign": str(experiment_remote),
            "epitope": spec["epitope"],
            "hotspots": epitope["hotspots"] if epitope else [56],
            "crop_clearance": clearances.get(key),
            "specification": spec,
            "control": epitope is None,
            "jobs": {},
            "glycans_modeled": False,
            "glycan_generation_guard": "asn_frame_excluded_volume_v1"
            if config.get("workflow") == "ig_followup_v1"
            else None,
        }
        (campaign / "manifest.json").write_text(json.dumps(record, indent=2) + "\n")
    record = {
        "campaign_id": args.campaign,
        "remote_project": str(project),
        "remote_campaign": str(remote),
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        **config,
        "jobs": {},
        "submission_intents": {},
        "bindcraft_revision": revision,
        "bindcraft_default_filter_sha256": file_hash(
            engine / "settings_filters/default_filters.json"
        ),
        "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "source_worktree_dirty": bool(subprocess.check_output(["git", "status", "--porcelain"])),
        "input_sha256": {
            str(p.relative_to(out)): file_hash(p) for p in out.rglob("*") if p.is_file()
        },
    }
    (out / "manifest.json").write_text(json.dumps(record, indent=2) + "\n")
    print(f"Prepared {out}: {len(config['experiments'])} configured experiments; no submission")


if __name__ == "__main__":
    main()
