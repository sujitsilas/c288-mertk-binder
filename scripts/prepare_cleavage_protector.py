#!/usr/bin/env python3
"""Freeze a bounded, explicitly exploratory MERTK stalk-protector campaign."""

import argparse
import hashlib
import json
import re
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from mertk.cleavage_protector import HOTSPOTS, prepare_targets
from mertk.protector_screen import disable_score_cutoffs
from mertk.provenance import file_hash
from mertk.target_policy import check_manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", required=True)
    parser.add_argument("--remote-project", required=True)
    args = parser.parse_args()
    if not re.fullmatch(r"round_[A-Za-z0-9_]+", args.campaign):
        raise ValueError("Invalid campaign identifier")
    out = Path("outputs/hpc") / args.campaign
    if out.exists():
        raise FileExistsError("Choose a fresh campaign; old inputs are immutable")
    config = json.loads(Path("configs/cleavage_protector_round.json").read_text())
    tools = json.loads(Path("configs/external_tools.json").read_text())
    engine = Path(tools["bindcraft"]["directory"])
    revision = subprocess.check_output(
        ["git", "-C", str(engine), "rev-parse", "HEAD"], text=True
    ).strip()
    if revision != tools["bindcraft"]["revision"]:
        raise ValueError("Unpinned BindCraft source")
    canonical = json.loads(Path("outputs/species_review_20261111/Q12866.json").read_text())[
        "sequence"
    ]["value"]
    source = Path("outputs/design/cleavage_protector_20261117")
    ref = out / "references"
    triage = prepare_targets(
        source / "AF_Q12866.pdb", Path("data/references/37NT/37NT.cif"), canonical, ref
    )
    for name in ["AF_Q12866.pdb", "AF_Q12866_pae.json", "alphafold_api.json", "fetch.json"]:
        shutil.copy2(source / name, ref / name)
    shutil.copy2("data/references/37NT/37NT.cif", ref / "37NT.cif")
    shutil.copy2("outputs/species_review_20261111/Q12866.json", ref / "canonical.json")
    workflow = out / "workflow"
    workflow.mkdir()
    for name in [
        "run_bounded_bindcraft.py",
        "check_glycan_runtime.py",
        "bindcraft_round.py",
        "check_course_outputs.py",
        "ig_followup.py",
        "protector_round.py",
        "run_protector_round.sh",
        "submit_protector_round.py",
        "submit_bindcraft_round.py",
    ]:
        shutil.copy2(Path("scripts/hpc") / name, workflow / name)
    shutil.copy2(__file__, workflow / "prepare_cleavage_protector.py")
    shutil.copytree("src/mertk", workflow / "mertk", ignore=shutil.ignore_patterns("__pycache__"))
    for name in [
        "external_tools.json",
        "active_target_policy.json",
        "cleavage_protector_round.json",
        "mertk_cleavage_protector.json",
    ]:
        shutil.copy2(Path("configs") / name, workflow / name)
    proof = Path(config["reused_positive_control"])
    if json.loads(proof.read_text())["passed"] is not True:
        raise ValueError("Positive-control evidence missing")
    shutil.copy2(proof, workflow / "positive_control_evidence.json")
    remote = args.remote_project.rstrip("/") + "/runs/" + args.campaign
    filters = (engine / "settings_filters/default_filters.json").read_bytes()
    if config.get("protector_screen") == "exploratory":
        settings = disable_score_cutoffs(json.loads(filters))
        filters = (json.dumps(settings, indent=2) + "\n").encode()
    for key, spec in config["experiments"].items():
        campaign = out / "campaigns" / key
        inp = campaign / "inputs"
        (inp / "bindcraft").mkdir(parents=True)
        shutil.copy2(ref / f"{spec['conformer']}_target.pdb", inp / "target.pdb")
        shutil.copy2(ref / f"{spec['conformer']}_constraints.json", inp / "glycan_constraints.json")
        settings = dict(
            chains="A",
            target_hotspot_residues=",".join(f"A{n}" for n in HOTSPOTS),
            design_path=f"{remote}/campaigns/{key}/bindcraft/designs",
            binder_name=key,
            starting_pdb=f"{remote}/campaigns/{key}/inputs/target.pdb",
            lengths=spec["lengths"],
            number_of_final_designs=401,
        )
        advanced = json.loads(
            (engine / "settings_advanced/default_4stage_multimer_hardtarget.json").read_text()
        )
        advanced.update(
            max_trajectories=spec["max_attempts"],
            num_seqs=20,
            max_mpnn_sequences=2,
            save_design_animations=False,
            save_design_trajectory_plots=True,
            save_mpnn_fasta=True,
            remove_unrelaxed_trajectory=False,
            remove_unrelaxed_complex=False,
            remove_binder_monomer=False,
            af_params_dir=args.remote_project + "/external/BindCraft",
        )
        for name, obj in [("target", settings), ("advanced", advanced)]:
            (inp / "bindcraft" / f"{name}.json").write_text(json.dumps(obj, indent=2) + "\n")
        (inp / "bindcraft/filters.json").write_bytes(filters)
        child = dict(
            campaign_id=args.campaign + "_" + key,
            remote_project=args.remote_project,
            remote_campaign=remote + "/campaigns/" + key,
            epitope="Stalk_P490S491",
            hotspots=HOTSPOTS,
            specification=spec,
            control=False,
            jobs={},
            glycan_generation_guard="asn_frame_excluded_volume_v1",
            glycans_modeled=False,
        )
        (campaign / "manifest.json").write_text(json.dumps(child, indent=2) + "\n")
    record = dict(
        campaign_id=args.campaign,
        remote_project=args.remote_project,
        remote_campaign=remote,
        created_at_utc=datetime.now(timezone.utc).isoformat(),
        **config,
        jobs={},
        submission_intents={},
        bindcraft_revision=revision,
        bindcraft_default_filter_sha256=file_hash(engine / "settings_filters/default_filters.json"),
        protector_filter_sha256=hashlib.sha256(filters).hexdigest(),
        target_triage=triage,
        canonical_target_sequence=canonical[92:526],
        source_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        input_sha256={str(p.relative_to(out)): file_hash(p) for p in out.rglob("*") if p.is_file()},
    )
    (out / "manifest.json").write_text(json.dumps(record, indent=2) + "\n")
    check_manifest(out / "manifest.json", generation=True)
    print(
        json.dumps(
            {
                "prepared": str(out),
                "frozen_files": len(record["input_sha256"]),
                "conformers": triage["conformers"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
