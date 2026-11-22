#!/usr/bin/env python3
"""Nominate a small, documented native-output cohort and prepare matched AF3 inputs."""

import argparse
import json
import re
import subprocess
from collections import defaultdict
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionException

from mertk.af3 import prepare_jobs
from mertk.candidates import import_complexes
from mertk.provenance import file_hash
from mertk.structural import rigid_fit, rmsd
from mertk.structures import read_design_chain


def chain_ids(path, sequence):
    parser = MMCIFParser(QUIET=True) if path.suffix == ".cif" else PDBParser(QUIET=True)
    structure = parser.get_structure("native", str(path))
    if len(structure) != 1 or len(structure[0]) != 2:
        raise ValueError("Not a single-model, two-chain complex")
    chains = {chain.id: read_design_chain(path, chain.id)[0] for chain in structure[0]}
    matches = [key for key, chain in chains.items() if chain.sequence == sequence]
    if len(matches) != 1:
        raise ValueError("No unique chain matching the full experimental target sequence")
    target = matches[0]
    binder = next(key for key in chains if key != target)
    return target, binder


def bindcraft_priority(path):
    parts = set(path.parts)
    priority = 0 if "Accepted" in parts else 1 if "Relaxed" in parts else 2
    return priority, str(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign", type=Path)
    parser.add_argument("--max-target-ca-rmsd", type=float, default=2.0)
    args = parser.parse_args()
    if not np.isfinite(args.max_target_ca_rmsd) or args.max_target_ca_rmsd <= 0:
        parser.error("--max-target-ca-rmsd must be finite and positive")
    campaign = args.campaign.resolve()
    manifest = json.loads((campaign / "manifest.json").read_text())
    target_path = campaign / "inputs/target.pdb"
    target, _ = read_design_chain(target_path, "A")
    output = campaign / "assessment"
    output.mkdir(exist_ok=False)
    policy = (
        "At most two complexes per generator; BindCraft accepted structures first, then "
        "relaxed trajectories, then other saved trajectory complexes, with one per actual "
        "trajectory seed. RFD4 generated complexes sorted by path. Require complete chains, "
        "exact target sequence and at least one inter-chain C-alpha contact within 8 A. "
        f"Target C-alpha RMSD after rigid alignment to experimental 37NT must be at most "
        f"{args.max_target_ca_rmsd:g} A. "
        "This runtime cohort is not an unbiased comparison of generator success rates."
    )
    selection = {"policy": policy, "selected": [], "rejected": []}
    files = {
        "bindcraft": sorted(
            [
                p
                for p in (campaign / "bindcraft/designs").rglob("*.pdb")
                if "Accepted" in p.parts or "Trajectory" in p.parts
            ],
            key=bindcraft_priority,
        ),
        "rfd4_proteina": sorted(
            (campaign / "rfd4_proteina/inference_samples").rglob("*_generated.cif")
        ),
    }
    rows, counts, seen, used_seeds = [], defaultdict(int), set(), set()
    for arm, paths in files.items():
        for path in paths:
            if counts[arm] >= 2:
                break
            digest = file_hash(path)
            if digest in seen:
                continue
            seen.add(digest)
            if arm == "bindcraft":
                match = re.search(r"_s(\d+)(?:_|$)", path.stem)
                if not match:
                    selection["rejected"].append({"path": str(path), "reason": "No native seed"})
                    continue
                seed = int(match[1])
                if seed in used_seeds:
                    continue
            else:
                seed = manifest["rfd4_seed"]
            run_record = output / f"{arm}_{counts[arm]}_run.json"
            try:
                target_id, binder_id = chain_ids(path, target.sequence)
                native_target, _ = read_design_chain(path, target_id)
                mobile, reference = native_target.backbone[:, 1], target.backbone[:, 1]
                rotation, translation = rigid_fit(mobile, reference)
                target_rmsd = rmsd(mobile @ rotation + translation, reference)
                if target_rmsd > args.max_target_ca_rmsd:
                    raise ValueError(
                        f"Target C-alpha RMSD {target_rmsd:.3f} A exceeds "
                        f"{args.max_target_ca_rmsd:g} A from experimental reference"
                    )
                run_record.write_text(
                    json.dumps(
                        {
                            "tool_revision": manifest["tools"][arm]["revision"],
                            "campaign_manifest_sha256": file_hash(campaign / "manifest.json"),
                            "generation_job": manifest["jobs"][f"{arm}_generate"],
                            "native_structure": str(path),
                            "selection_policy": policy,
                            "native_acceptance": "accepted"
                            if "Accepted" in path.parts
                            else "not established",
                        },
                        indent=2,
                    )
                    + "\n"
                )
                imported = import_complexes(
                    [path],
                    output / f"{arm}_{counts[arm]}",
                    arm=arm,
                    campaign_id=manifest["campaign_id"],
                    epitope_id="unreviewed_surface",
                    target_reference=target_path,
                    run_record=run_record,
                    seed=seed,
                    binder_chain=binder_id,
                    target_chain=target_id,
                )
            except (ValueError, KeyError, PDBConstructionException) as error:
                selection["rejected"].append({"path": str(path), "reason": str(error)})
                continue
            rows.extend(imported)
            counts[arm] += 1
            if arm == "bindcraft":
                used_seeds.add(seed)
            selection["selected"].append(
                {
                    "arm": arm,
                    "path": str(path),
                    "seed": seed,
                    "candidate_id": imported[0]["candidate_id"],
                    "target_ca_rmsd": target_rmsd,
                }
            )
    selection["counts"] = dict(counts)
    (output / "selection.json").write_text(json.dumps(selection, indent=2) + "\n")
    candidates = output / "candidates.jsonl"
    candidates.write_text("".join(json.dumps(row) + "\n" for row in rows))
    if any(counts[arm] == 0 for arm in files):
        raise SystemExit("AF3 cohort not prepared: at least one generator has no usable complex.")
    protocol = json.loads(Path("configs/af3_pilot_single_sequence.json").read_text())
    plan = prepare_jobs(candidates, protocol, output / "af3")
    print(f"Prepared {len(plan['jobs'])} AF3 candidates; {dict(counts)}")
    result = subprocess.run(
        [
            "sacct",
            "--allocations",
            "-nP",
            "-j",
            ",".join(job["job_id"] for job in manifest["jobs"].values()),
            "-o",
            "JobID,State,ExitCode,Elapsed",
        ],
        text=True,
        capture_output=True,
    )
    (output / "generation_accounting.txt").write_text(result.stdout + result.stderr)


if __name__ == "__main__":
    main()
