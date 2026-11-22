#!/usr/bin/env python3
"""Prepare a fresh, bounded Hoffman2 pilot from the verified experimental reference."""

import argparse
import json
import re
import shutil
import subprocess
from pathlib import Path, PurePosixPath

from mertk.provenance import file_hash
from mertk.structures import read_design_chain


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", required=True)
    parser.add_argument("--remote-project", required=True)
    parser.add_argument("--reference", type=Path, default=Path("data/references/37NT"))
    parser.add_argument("--output-root", type=Path, default=Path("outputs/hpc"))
    args = parser.parse_args()
    parser.error(
        "Unrestricted MERTK surface exploration is retired. "
        "Use prepare_bindcraft_round.py with the active Ig1/Ig2 targets."
    )
    if not re.fullmatch(r"[A-Za-z0-9_-]+", args.campaign):
        parser.error("Campaign must contain only letters, digits, underscores and hyphens")
    project = PurePosixPath(args.remote_project)
    if not project.is_absolute():
        parser.error("Remote project must be an absolute path")
    source = json.loads((args.reference / "provenance.json").read_text())
    target, mapping = read_design_chain(args.reference / "target.pdb", "A")
    if source["pdb_id"] != "37NT" or source["author_chain"] != "A":
        raise ValueError("This pilot is defined for experimental 37NT chain A")
    if [r["author_number"] for r in mapping] != list(range(93, 485)):
        raise ValueError("Pilot requires the complete observed A93-484 construct")
    if file_hash(args.reference / "target.pdb") != source["target_sha256"]:
        raise ValueError("Exported target changed since reference preparation")
    if file_hash(args.reference / "37NT.cif") != source["source_sha256"]:
        raise ValueError("Original experimental mmCIF changed since reference preparation")
    if target.sequence != source["sequence"]:
        raise ValueError("Target sequence differs from the reference provenance")
    output = args.output_root / args.campaign
    if output.exists():
        raise FileExistsError("Choose a fresh campaign; existing runs are never overwritten")
    remote = project / "runs" / args.campaign
    tools = json.loads(Path("configs/external_tools.json").read_text())
    bindcraft = Path(tools["bindcraft"]["directory"])
    advanced = json.loads(
        (bindcraft / "settings_advanced/default_4stage_multimer.json").read_text()
    )
    advanced.update(
        max_trajectories=2,
        save_design_animations=False,
        save_design_trajectory_plots=True,
        save_mpnn_fasta=True,
        remove_unrelaxed_trajectory=False,
        remove_unrelaxed_complex=False,
        remove_binder_monomer=False,
        af_params_dir=str(project / "external/BindCraft"),
    )
    settings = {
        "design_path": str(remote / "bindcraft/designs"),
        "binder_name": "MERTK_surface_pilot",
        "starting_pdb": str(remote / "inputs/target.pdb"),
        "chains": "A",
        "target_hotspot_residues": None,
        "lengths": [80, 80],
        "number_of_final_designs": 1,
    }
    inputs = output / "inputs"
    (inputs / "bindcraft").mkdir(parents=True)
    for name in ("target.pdb", "37NT.cif"):
        shutil.copy2(args.reference / name, inputs / name)
    shutil.copy2(args.reference / "provenance.json", inputs / "reference_provenance.json")
    shutil.copy2(
        bindcraft / "settings_filters/default_filters.json", inputs / "bindcraft/filters.json"
    )
    for name, data in (("target", settings), ("advanced", advanced)):
        (inputs / f"bindcraft/{name}.json").write_text(json.dumps(data, indent=2) + "\n")
    (inputs / "rfd4.yaml").write_text(
        "# Surface exploration only; no selected epitope or agonism claim.\n"
        "mertk_surface_pilot:\n  input: ./37NT.cif\n  input_as_annotated_cif: false\n"
        '  contig: "80, /0, A93-484"\n  conditions:\n    C_CRDGRP:\n      1:\n'
        "        - select: \"sel('A/*/*/[N,CA,C,O]') and not has_nan_coord()\"\n"
        "    C_SEQ:\n      true:\n        - select: \"sel('A')\"\n"
    )
    manifest = {
        "campaign_id": args.campaign,
        "remote_project": str(project),
        "remote_campaign": str(remote),
        "purpose": "Runtime and unreviewed surface exploration",
        "reference": "37NT chain A residues 93-484",
        "binder_length": 80,
        "epitope_selected": False,
        "glycans_modeled": False,
        "bindcraft_max_trajectories": 2,
        "bindcraft_seed": "Recover actual seeds from native logs",
        "rfd4_replicates": 2,
        "rfd4_seed": 101,
        "rfd4_recycle_strategy": "medium/front_loaded",
        "af3_status": "Not run",
        "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "tools": tools,
        "jobs": {},
        "input_sha256": {
            str(p.relative_to(output)): file_hash(p) for p in inputs.rglob("*") if p.is_file()
        },
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Prepared {output}; no jobs submitted")


if __name__ == "__main__":
    main()
