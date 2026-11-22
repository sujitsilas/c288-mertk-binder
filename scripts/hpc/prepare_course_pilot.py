#!/usr/bin/env python3
"""Prepare a fresh, small two-generator run plus a gated upstream RFD4 control."""

import argparse
import json
import re
import shutil
import subprocess
from pathlib import Path, PurePosixPath

import numpy as np
from Bio.PDB import PDBIO, PDBParser, Select

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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", required=True)
    parser.add_argument("--remote-project", required=True)
    epitopes = json.loads(Path("configs/course_epitopes.json").read_text())
    parser.add_argument("--epitope", choices=sorted(epitopes), default="Ig1_E152")
    parser.add_argument("--output-root", type=Path, default=Path("outputs/hpc"))
    # Sampling depth is the dominant determinant of whether any design is accepted.
    # Upstream guidance: ~100 trajectories for an easy target, 1,000-10,000 for a hard
    # one. The November 12 pilot ran 4 and 6 attempts and accepted nothing, which is the
    # expected outcome at that depth and says nothing about the epitope. Do not lower
    # these to economise: cutting num_seqs starves the only stage that converts a good
    # backbone into a sequence that survives AF2 revalidation.
    parser.add_argument(
        "--max-attempts",
        type=int,
        default=200,
        help="Admitted design-loop iterations, including rejected trials",
    )
    parser.add_argument(
        "--num-seqs", type=int, default=20, help="MPNN sequences sampled per completed trajectory"
    )
    parser.add_argument(
        "--max-mpnn-sequences",
        type=int,
        default=2,
        help="Accepted MPNN variants kept per trajectory",
    )
    parser.add_argument("--final-designs", type=int, default=4)
    parser.add_argument(
        "--length-range",
        type=int,
        nargs=2,
        default=[65, 130],
        metavar=("MIN", "MAX"),
        help="Binder lengths are sampled; a single fixed length removes "
        "a free and effective axis of diversity",
    )
    args = parser.parse_args()
    if args.num_seqs < 8:
        parser.error(
            "num_seqs below 8 starves MPNN resampling; the pilot's 40 sampled sequences "
            "all failed AF2 revalidation, which argues for more sampling, not less"
        )
    if args.length_range[0] > args.length_range[1]:
        parser.error("Length range must be ascending")
    epitope = epitopes[args.epitope]
    require_ig_target(args.epitope, epitope["hotspots"])
    start, end = epitope["crop"]
    hotspots = epitope["hotspots"]
    if not re.fullmatch(r"[A-Za-z0-9_-]+", args.campaign):
        parser.error("Invalid campaign identifier")
    project = PurePosixPath(args.remote_project)
    if not project.is_absolute():
        parser.error("Remote project must be absolute")
    out = args.output_root / args.campaign
    if out.exists():
        raise FileExistsError("Existing campaign is immutable; choose a fresh identifier")
    reference = Path("data/references/37NT")
    provenance = json.loads((reference / "provenance.json").read_text())
    for name, key in [("37NT.cif", "source_sha256"), ("target.pdb", "target_sha256")]:
        if file_hash(reference / name) != provenance[key]:
            raise ValueError("Experimental reference checksum changed")
    tools = json.loads(Path("configs/external_tools.json").read_text())
    bindcraft = Path(tools["bindcraft"]["directory"])
    rfd4 = Path(tools["rfd4_proteina"]["directory"])
    for name, folder in [("bindcraft", bindcraft), ("rfd4_proteina", rfd4)]:
        actual = subprocess.check_output(
            ["git", "-C", str(folder), "rev-parse", "HEAD"], text=True
        ).strip()
        if actual != tools[name]["revision"]:
            raise ValueError(f"Unpinned {name} checkout")
    inputs = out / "inputs"
    (inputs / "bindcraft").mkdir(parents=True)
    for name in ["37NT.cif", "provenance.json"]:
        shutil.copy2(reference / name, inputs / name)
    shutil.copy2(reference / "target.pdb", inputs / "full_target.pdb")
    writer = PDBIO()
    writer.set_structure(PDBParser(QUIET=True).get_structure("native", reference / "target.pdb"))
    writer.save(str(inputs / "target.pdb"), select=Crop(start, end))
    target, mapping = read_design_chain(inputs / "target.pdb", "A")
    full, _ = read_design_chain(reference / "target.pdb", "A")
    assert [r["author_number"] for r in mapping] == list(range(start, end + 1))
    assert target.sequence == full.sequence[start - 93 : end - 93 + 1]
    assert np.allclose(target.backbone, full.backbone[start - 93 : end - 93 + 1])
    (inputs / "residue_map.json").write_text(json.dumps(mapping, indent=2) + "\n")
    native = PDBParser(QUIET=True).get_structure("crop", inputs / "target.pdb")[0]["A"]
    native_reference = {
        "source": "37NT; original author residue numbers; coordinates rounded by PDB export",
        "residues": [
            {
                "author_number": r.id[1],
                "residue_name": r.resname,
                "backbone": {name: r[name].coord.tolist() for name in ["N", "CA", "C", "O"]},
            }
            for r in native
        ],
        "hotspot_atoms": [
            [n, a.name]
            for n in hotspots
            for a in native[n]
            if a.element != "H"
            and a.name not in {"N", "CA", "C", "O"}
            and np.isfinite(a.coord).all()
        ],
    }
    (inputs / "target_reference.json").write_text(json.dumps(native_reference, indent=2) + "\n")
    # A design crop must carry the epitope's whole structural neighbourhood. The November
    # 12 Patch A crop began at A286 and so removed I283 and P284, 3.2-3.7 A from the
    # epitope centroid; binders then expanded into that void and collided with Ig2 R210
    # once the intact receptor was restored for nomination. Designing and judging on
    # different molecules guarantees rejections that look like epitope failure.
    intact = PDBParser(QUIET=True).get_structure("intact", reference / "target.pdb")[0]["A"]
    intact_residues = {r.id[1]: r for r in intact if r.id[0] == " "}
    centroid = np.mean(
        [intact_residues[n]["CA"].coord for n in hotspots if "CA" in intact_residues[n]], axis=0
    )
    severed = sorted(
        (float(min(np.linalg.norm(a.coord - centroid) for a in r if a.element != "H")), n)
        for n, r in intact_residues.items()
        if not (start <= n <= end)
    )
    intruding = [(round(d, 2), n) for d, n in severed if d < 20.0]
    if intruding:
        raise ValueError(
            f"Crop A{start}-A{end} excludes {len(intruding)} residues within 20 A of the "
            f"{args.epitope} centroid, closest {intruding[:5]}; widen the crop to a domain "
            "boundary that contains the epitope neighbourhood"
        )
    native_reference["crop_clearance"] = {
        "nearest_excluded_residue_angstroms": round(severed[0][0], 2) if severed else None,
        "policy": "no excluded residue within 20 A of the epitope centroid",
    }
    (inputs / "target_reference.json").write_text(json.dumps(native_reference, indent=2) + "\n")
    remote = project / "runs" / args.campaign
    advanced = json.loads(
        (bindcraft / "settings_advanced/default_4stage_multimer.json").read_text()
    )
    advanced.update(
        max_trajectories=args.max_attempts,
        num_seqs=args.num_seqs,
        max_mpnn_sequences=args.max_mpnn_sequences,
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
        "binder_name": f"MERTK_{args.epitope}_pilot",
        "starting_pdb": str(remote / "inputs/target.pdb"),
        "chains": "A",
        "target_hotspot_residues": ",".join(f"A{n}" for n in hotspots),
        "lengths": list(args.length_range),
        "number_of_final_designs": args.final_designs,
    }
    for name, data in [("advanced", advanced), ("target", settings)]:
        (inputs / f"bindcraft/{name}.json").write_text(json.dumps(data, indent=2) + "\n")
    shutil.copy2(
        bindcraft / "settings_filters/default_filters.json", inputs / "bindcraft/filters.json"
    )
    (inputs / "rfd4.yaml").write_text(
        f"# Exploratory {args.epitope}; protein-only generation, intact glycan review required.\n"
        "mertk_design:\n  input: ./37NT.cif\n  input_as_annotated_cif: false\n"
        f'  contig: "80, /0, A{start}-{end}"\n  conditions:\n    C_CRDGRP:\n      1:\n'
        "        - select: \"sel('A/*/*/[N,CA,C,O]') and not has_nan_coord()\"\n"
        "    C_SEQ:\n      true:\n        - select: \"sel('A')\"\n"
        "    C_HOT:\n      true:\n"
        f"        - select: \"sel('A/*/[{','.join(map(str, hotspots))}]') and not sel('A/*/*/[N,CA,C,O]') and not has_nan_coord()\"\n"
    )
    example = rfd4 / "docs/examples/inputs/5o45_cropped.cif.zst"
    raw = subprocess.check_output(["zstd", "-dc", str(example)])
    (inputs / "pdl1_control.cif").write_bytes(raw)
    control = (rfd4 / "docs/how-to/input_prep/examples/condition_spec_ppi.yaml").read_text()
    control = control.replace("../../../examples/inputs/5o45_cropped.cif.zst", "./pdl1_control.cif")
    (inputs / "rfd4_control.yaml").write_text(control)
    manifest = {
        "campaign_id": args.campaign,
        "remote_project": str(project),
        "remote_campaign": str(remote),
        "purpose": "Bounded epitope feasibility pilot; RFD4 control gates its MERTK generation; no efficacy or tool-superiority claim",
        "epitope": args.epitope,
        "epitope_description": epitope["description"],
        "reference": f"37NT A{start}-{end} design crop; intact A93-484 for context review",
        "binder_length_range": list(args.length_range),
        "hotspots": hotspots,
        "target_residues": len(mapping),
        "epitope_selected": "provisional exploratory anchors",
        "glycans_modeled": False,
        "bindcraft_max_trajectories": args.max_attempts,
        "bindcraft_max_attempts": args.max_attempts,
        "bindcraft_num_sampled_sequences_per_completed_trajectory": args.num_seqs,
        "bindcraft_attempt_policy": "Project launcher counts every admitted design-loop iteration, including upstream-rejected trials; max_trajectories alone only counts relaxed structures",
        "rfd4_replicates": 2,
        "rfd4_seed": 202,
        "rfd4_control_seed": 101,
        "rfd4_recycle_strategy": "medium/front_loaded",
        "af3_status": "Not submitted; requires real candidates passing geometry and native-context review",
        "authorized_account": "pi_pallard",
        "partition": "campus",
        "qos": "campus24",
        "maximum_allocated_gpu_hours_this_tranche": 7,
        "budgets": {
            "bindcraft_generate": "04:00:00",
            "rfd4_control": "01:00:00",
            "rfd4_proteina_generate": "02:00:00",
        },
        "tools": tools,
        "jobs": {},
        "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "input_sha256": {
            str(p.relative_to(out)): file_hash(p) for p in inputs.rglob("*") if p.is_file()
        },
        "control_policy": "One upstream PD-L1 PPI example, unchanged except input path; verify target sequence, <=2 A target CA RMSD and non-catastrophic adjacent CA geometry before starting MERTK",
        "quality_gates": {
            "max_target_ca_rmsd_A": 2.0,
            "adjacent_ca_range_A": [2.7, 4.5],
            "max_intact_protein_clash_distance_A": 2.0,
            "min_observed_glycan_clearance_A": 2.5,
        },
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(
        f"Prepared {out}; {len(mapping)}-residue native target, provisional {args.epitope}, two bounded arms; no submission"
    )


if __name__ == "__main__":
    main()
