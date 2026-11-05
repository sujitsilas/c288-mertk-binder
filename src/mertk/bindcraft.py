"""Prepare BindCraft jobs and import completed complexes as application holdouts."""

import hashlib
import json
import re
import shlex
import shutil
import subprocess
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

from mertk.data import AA, BACKBONE, ELEMENTS, Chain, prepare_pair
from mertk.training import file_hash


def read_design_chain(path, chain_id):
    model = PDBParser(QUIET=True).get_structure("design", str(path))[0]
    if chain_id not in model:
        raise ValueError(f"Chain {chain_id!r} is absent from {path}")
    residues = [r for r in model[chain_id] if r.id[0] == " "]
    if not residues:
        raise ValueError("Selected chain has no protein residues")
    sequence, backbone, oxygen, atoms, features, mapping = [], [], [], [], [], []
    previous = None
    for residue in residues:
        aa = seq1(residue.resname, undef_code="X")
        if aa not in AA or not all(name in residue for name in ["N", "CA", "C", "O"]):
            raise ValueError(
                "Application PDB chains must have standard residues and complete N/CA/C/O backbones"
            )
        number = residue.id[1]
        if previous is not None and number not in {previous, previous + 1}:
            raise ValueError(
                "PDB chain has numbering gaps; prepare a continuous, explicitly renumbered construct first"
            )
        previous = number
        sequence.append(aa)
        backbone.append([residue[name].coord for name in ["N", "CA", "C"]])
        oxygen.append(residue["O"].coord)
        mapping.append(
            {
                "position": len(sequence),
                "author_number": number,
                "insertion_code": residue.id[2].strip(),
            }
        )
        for atom in residue:
            if atom.element in {"H", "D"} or (atom.occupancy or 0) <= 0:
                continue
            if not np.isfinite(atom.coord).all():
                raise ValueError("Non-finite PDB coordinates")
            atoms.append(atom.coord)
            features.append(
                [ELEMENTS.get(atom.element, 5), AA.index(aa), int(atom.name in BACKBONE)]
            )
    backbone = np.array(backbone, dtype=np.float32)
    oxygen = np.array(oxygen, dtype=np.float32)
    if not np.isfinite(backbone).all() or not np.isfinite(oxygen).all():
        raise ValueError("Non-finite backbone")
    return Chain(
        chain_id,
        "".join(sequence),
        backbone,
        np.array(atoms, dtype=np.float32),
        np.array(features, dtype=np.int64),
        oxygen,
    ), mapping


def validate_hotspots(hotspots, chain_id, mapping):
    if hotspots is None:
        return None
    positions = {r["author_number"] for r in mapping}
    for token in hotspots.split(","):
        match = re.fullmatch(r"([A-Za-z]?)(\d+)(?:-(\d+))?", token.strip())
        if not match or (match[1] and match[1] != chain_id):
            raise ValueError(
                "Hotspots must be residue numbers/ranges on the selected chain, e.g. A25,A30-33"
            )
        start, end = int(match[2]), int(match[3] or match[2])
        if start > end or not set(range(start, end + 1)) <= positions:
            raise ValueError(f"Hotspot residue is absent: {token}")
    return hotspots


def prepare_job(
    target_pdb,
    target_chain,
    output,
    lengths=(65, 120),
    designs=10,
    hotspots=None,
    name="MerTK",
    max_length=512,
):
    output = Path(output).resolve()
    target_pdb = Path(target_pdb).resolve()
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"Choose an empty job directory: {output}")
    if not re.fullmatch(r"[A-Za-z0-9_-]+", name):
        raise ValueError("Use an alphanumeric binder name")
    if not 1 <= lengths[0] <= lengths[1] or designs < 1:
        raise ValueError("Invalid binder length range or design count")
    target, mapping = read_design_chain(target_pdb, target_chain)
    if len(target.sequence) + lengths[1] + 3 > max_length:
        raise ValueError(
            "Target plus maximum binder length exceeds the redesign token budget; choose a smaller documented construct"
        )
    hotspots = validate_hotspots(hotspots, target_chain, mapping)
    spec = json.loads(Path("configs/external_tools.json").read_text())["bindcraft"]
    root = Path(spec["directory"]).resolve()
    if not root.exists():
        raise FileNotFoundError("Run scripts/setup_external.py --tool bindcraft")
    revision = subprocess.check_output(
        ["git", "-C", str(root), "rev-parse", "HEAD"], text=True
    ).strip()
    dirty = subprocess.check_output(
        ["git", "-C", str(root), "status", "--porcelain", "--untracked-files=no"], text=True
    )
    if revision != spec["revision"] or dirty:
        raise ValueError("BindCraft checkout differs from the pinned revision")
    output.mkdir(parents=True, exist_ok=True)
    shutil.copy2(target_pdb, output / "target.pdb")
    filters = root / "settings_filters/default_filters.json"
    advanced = root / "settings_advanced/default_4stage_multimer.json"
    shutil.copy2(filters, output / "filters.json")
    shutil.copy2(advanced, output / "advanced.json")
    settings = {
        "design_path": str(output / "designs"),
        "binder_name": name,
        "starting_pdb": str(output / "target.pdb"),
        "chains": target_chain,
        "target_hotspot_residues": hotspots,
        "lengths": list(lengths),
        "number_of_final_designs": designs,
    }
    (output / "target.json").write_text(json.dumps(settings, indent=2) + "\n")
    command = [
        "conda",
        "run",
        "--no-capture-output",
        "-n",
        "BindCraft",
        "python",
        "-u",
        str(root / "bindcraft.py"),
        "--settings",
        str(output / "target.json"),
        "--filters",
        str(output / "filters.json"),
        "--advanced",
        str(output / "advanced.json"),
    ]
    (output / "run.sh").write_text(
        "#!/usr/bin/env bash\nset -euo pipefail\ncd "
        + shlex.quote(str(root))
        + "\nexec "
        + shlex.join(command)
        + "\n"
    )
    provenance = {
        "tool": spec,
        "target_sha256": file_hash(output / "target.pdb"),
        "target_chain": target_chain,
        "target_sequence": target.sequence,
        "residue_mapping": mapping,
        "settings": settings,
        "filters_sha256": file_hash(output / "filters.json"),
        "advanced_sha256": file_hash(output / "advanced.json"),
        "launched": False,
    }
    (output / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    return output


def import_complexes(paths, output, binder_chain="B", target_chain="A", max_length=512):
    if binder_chain == target_chain:
        raise ValueError("Binder and target must be different chains")
    output = Path(output).resolve()
    if (output / "examples.jsonl").exists():
        raise FileExistsError("Application manifest already exists; choose another output")
    output.mkdir(parents=True, exist_ok=True)
    rows, provenance, fasta = [], [], []
    for source in paths:
        source = Path(source).resolve()
        binder, binder_mapping = read_design_chain(source, binder_chain)
        target, target_mapping = read_design_chain(source, target_chain)
        if len(binder.sequence) + len(target.sequence) + 3 > max_length:
            raise ValueError(f"Complex exceeds the token budget: {source}")
        example = prepare_pair(binder, target)
        if not example["design_mask"].any():
            raise ValueError(f"No binder interface contacts: {source}")
        source_hash = file_hash(source)
        identifier = re.sub(r"[^A-Za-z0-9_-]", "_", source.stem) + "_" + source_hash[:10]
        if any(row["id"] == identifier for row in rows):
            raise ValueError("Duplicate input complex")
        np.savez_compressed(output / f"{identifier}.npz", **example)
        rows.append(
            {
                "id": identifier,
                "path": f"{identifier}.npz",
                "split": "application_holdout",
                "component": "application_target_"
                + hashlib.sha256(target.sequence.encode()).hexdigest()[:20],
                "pdb_id": "application_" + source_hash[:10],
                "binder_key": "seq_" + hashlib.sha256(binder.sequence.encode()).hexdigest()[:20],
                "target_key": "seq_" + hashlib.sha256(target.sequence.encode()).hexdigest()[:20],
                "source_pdb": str(source),
                "source_sha256": source_hash,
                "binder_chain": binder_chain,
                "target_chain": target_chain,
            }
        )
        provenance.append(
            {
                "id": identifier,
                "source": str(source),
                "source_sha256": source_hash,
                "binder_chain": binder_chain,
                "target_chain": target_chain,
                "binder_mapping": binder_mapping,
                "target_mapping": target_mapping,
                "baseline_sequence": binder.sequence,
            }
        )
        fasta.append(f">{identifier} source=BindCraft_original\n{binder.sequence}\n")
    if not rows:
        raise ValueError("No input complexes")
    (output / "examples.jsonl").write_text("".join(json.dumps(row) + "\n" for row in rows))
    (output / "bindcraft_original.fasta").write_text("".join(fasta))
    (output / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    return rows


def export_candidates(rows, results, output):
    """Keep each original control alongside redesigns for a shared external scorer."""
    candidates, contexts = [], {}
    for row in rows:
        if file_hash(row["source_pdb"]) != row["source_sha256"]:
            raise ValueError("Original BindCraft complex changed after import")
        with np.load(row["path"], allow_pickle=False) as example:
            context = {
                "complex_id": row["id"],
                "source_pdb": row["source_pdb"],
                "source_sha256": row["source_sha256"],
                "binder_chain": row["binder_chain"],
                "target_chain": row["target_chain"],
                "target_sequence": str(example["target_sequence"]),
                "design_positions_1based": (np.flatnonzero(example["design_mask"]) + 1).tolist(),
                "structure_evaluation": "pending",
            }
            candidates.append(
                {
                    **context,
                    "candidate_id": row["id"] + ":bindcraft_original",
                    "arm": "bindcraft_original",
                    "seed": None,
                    "binder_sequence": str(example["binder_sequence"]),
                }
            )
            contexts[row["id"]] = context
    for result in results:
        if result.get("reference_kind") != "bindcraft_seed":
            raise ValueError("Expected an application comparison with BindCraft seed references")
        for case in result["cases"]:
            for index, sequence in enumerate(case["sequences"]):
                candidates.append(
                    {
                        **contexts[case["id"]],
                        "candidate_id": f"{case['id']}:{result['arm']}:{result['seed']}:{index}",
                        "arm": result["arm"],
                        "seed": result["seed"],
                        "binder_sequence": sequence,
                    }
                )
    Path(output).write_text("".join(json.dumps(candidate) + "\n" for candidate in candidates))
    return candidates
