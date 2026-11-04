"""Prepare fixed-backbone chain-pair examples without binder side-chain features."""

import gzip
import hashlib
import io
import json
import re
from dataclasses import dataclass
from itertools import permutations
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SeqUtils import seq1

AA = "ACDEFGHIKLMNPQRSTVWY"
ELEMENTS = {"C": 0, "N": 1, "O": 2, "S": 3, "P": 4}
BACKBONE = {"N", "CA", "C", "O", "OXT"}


@dataclass
class Chain:
    name: str
    sequence: str
    backbone: np.ndarray
    atoms: np.ndarray
    atom_features: np.ndarray


def read_chains(path):
    path = Path(path)
    text = gzip.open(path, "rt").read() if path.suffix == ".gz" else path.read_text()
    metadata = MMCIF2Dict(io.StringIO(text))
    sequences = dict(
        zip(
            metadata["_entity_poly.entity_id"],
            metadata["_entity_poly.pdbx_seq_one_letter_code_can"],
        )
    )
    entities = dict(zip(metadata["_struct_asym.id"], metadata["_struct_asym.entity_id"]))
    model = MMCIFParser(QUIET=True, auth_chains=False, auth_residues=False).get_structure(
        path.stem, io.StringIO(text)
    )[0]
    chains = []
    for chain in model:
        # Expanded assembly label_asym_ids must be used, never ambiguous author IDs.
        entity = entities.get(chain.id)
        if entity not in sequences:
            continue
        sequence = re.sub(r"\s", "", sequences[entity])
        sequence = re.sub(r"\([^)]*\)", "X", sequence)
        sequence = "".join(a if a in AA else "X" for a in sequence)
        coords = np.full((len(sequence), 3, 3), np.nan, dtype=np.float32)
        atoms, features = [], []
        for residue in chain:
            index = residue.id[1] - 1  # label_seq_id indexes entity_poly sequence
            aa = seq1(residue.resname, undef_code="X")
            if not 0 <= index < len(sequence) or aa not in AA or sequence[index] != aa:
                continue
            if all(name in residue for name in ["N", "CA", "C"]):
                coords[index] = [residue[name].coord for name in ["N", "CA", "C"]]
            for atom in residue:
                if atom.element in {"H", "D"} or (atom.occupancy or 0) <= 0:
                    continue
                if not np.isfinite(atom.coord).all():
                    continue
                atoms.append(atom.coord)
                features.append(
                    [ELEMENTS.get(atom.element, 5), AA.index(aa), int(atom.name in BACKBONE)]
                )
        if atoms and np.isfinite(coords).all(axis=(1, 2)).sum() >= 16:
            chains.append(
                Chain(
                    chain.id,
                    sequence,
                    coords,
                    np.array(atoms, dtype=np.float32),
                    np.array(features, dtype=np.int64),
                )
            )
    return chains


def neighbors(backbone, atoms, atom_features, k=32, radius=12.0):
    """Only target atoms enter this function. Distances are from binder CA atoms."""
    if k < 1 or radius <= 0:
        raise ValueError("k and radius must be positive")
    valid = np.isfinite(backbone).all(axis=(1, 2))
    delta = np.nan_to_num(backbone[:, 1, None, :]) - atoms[None, :, :]
    distances = np.linalg.norm(delta, axis=-1)
    count = min(k, len(atoms))
    selected = np.argsort(distances, axis=-1, kind="stable")[:, :count]
    selected_distances = np.take_along_axis(distances, selected, axis=1)
    mask = (selected_distances <= radius) & valid[:, None]
    xyz = np.zeros((len(backbone), k, 3), dtype=np.float32)
    features = np.zeros((len(backbone), k, 3), dtype=np.int64)
    padded_mask = np.zeros((len(backbone), k), dtype=bool)
    xyz[:, :count], features[:, :count], padded_mask[:, :count] = (
        atoms[selected],
        atom_features[selected],
        mask,
    )
    return xyz, features, padded_mask, np.where(valid, distances.min(axis=1), np.inf)


def prepare_pair(binder, target, k=32, radius=12.0, interface_cutoff=8.0):
    if not 0 < interface_cutoff <= radius:
        raise ValueError("interface_cutoff must be positive and <= neighbor radius")
    xyz, features, neighbor_mask, distance = neighbors(
        binder.backbone, target.atoms, target.atom_features, k, radius
    )
    # Interface membership uses binder backbone only, independent of its native residue identity.
    interface = (distance <= interface_cutoff) & np.array([a in AA for a in binder.sequence])
    return {
        "binder_sequence": np.array(binder.sequence),
        "target_sequence": np.array(target.sequence),
        "binder_backbone": binder.backbone,
        "target_backbone": target.backbone,
        "neighbor_xyz": xyz,
        "neighbor_features": features,
        "neighbor_mask": neighbor_mask,
        "design_mask": interface,
    }


def build_dataset(raw_manifest, output, max_length=512, min_interface=5, k=32, radius=12.0):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    source = json.loads(Path(raw_manifest).read_text())
    records, rejected, sequences = [], [], {}
    for item in source["assemblies"]:
        chains = read_chains(item["path"])
        for binder, target in permutations(chains, 2):
            # Homomers give away the redesign labels through the visible target sequence.
            if binder.sequence == target.sequence:
                continue
            if len(binder.sequence) + len(target.sequence) + 3 > max_length:
                continue
            example = prepare_pair(binder, target, k=k, radius=radius)
            if example["design_mask"].sum() < min_interface:
                continue
            identifier = f"{item['assembly']}_{binder.name}_{target.name}"
            path = output / f"{identifier}.npz"
            np.savez_compressed(path, **example)
            keys = []
            for chain in (binder, target):
                key = "seq_" + hashlib.sha256(chain.sequence.encode()).hexdigest()[:20]
                sequences[key] = chain.sequence
                keys.append(key)
            records.append(
                {
                    "id": identifier,
                    "path": path.name,
                    "assembly": item["assembly"],
                    "pdb_id": item["assembly"].split("-")[0],
                    "binder_chain": binder.name,
                    "target_chain": target.name,
                    "binder_key": keys[0],
                    "target_key": keys[1],
                    "split": "unsplit",
                    "interface_residues": int(example["design_mask"].sum()),
                }
            )
        if not any(r["assembly"] == item["assembly"] for r in records):
            rejected.append(item["assembly"])
    (output / "examples.jsonl").write_text("".join(json.dumps(r) + "\n" for r in records))
    fasta = "".join(f">{key}\n{sequence}\n" for key, sequence in sorted(sequences.items()))
    fasta += Path(source["tam_seeds"]).read_text()
    (output / "chains.fasta").write_text(fasta)
    (output / "preparation.json").write_text(
        json.dumps(
            {
                "source_manifest": str(Path(raw_manifest).resolve()),
                "max_length": max_length,
                "min_interface": min_interface,
                "k": k,
                "radius": radius,
                "interface_cutoff": 8.0,
                "examples": len(records),
                "no_examples": rejected,
            },
            indent=2,
        )
        + "\n"
    )
    return records
