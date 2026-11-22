"""Import externally designed complexes for a common structural assessment."""

import json
from pathlib import Path

import numpy as np

from mertk.provenance import file_hash
from mertk.structures import read_design_chain


def import_complexes(
    paths,
    output,
    *,
    arm,
    campaign_id,
    epitope_id,
    target_reference,
    run_record,
    seed,
    binder_chain="B",
    target_chain="A",
):
    if binder_chain == target_chain:
        raise ValueError("Binder and target must be different chains")
    if not all(isinstance(s, str) and s.strip() for s in (arm, campaign_id, epitope_id)):
        raise ValueError("Provide nonempty arm, campaign and epitope identifiers")
    if type(seed) is not int or seed < 0:
        raise ValueError("Seed must be a nonnegative integer")
    output = Path(output).resolve()
    if output.exists() and any(output.iterdir()):
        raise FileExistsError("Choose an empty import directory")
    paths = [Path(p).resolve() for p in paths]
    if not paths:
        raise ValueError("No input complexes")
    run = json.loads(Path(run_record).read_text())
    if not isinstance(run, dict) or not run.get("tool_revision"):
        raise ValueError("Run record must identify the tool_revision")
    target_hash, run_hash = file_hash(target_reference), file_hash(run_record)
    prepared, seen = [], set()
    for index, path in enumerate(paths):
        source_hash = file_hash(path)
        if source_hash in seen:
            raise ValueError("Duplicate input complex")
        seen.add(source_hash)
        binder, binder_mapping = read_design_chain(path, binder_chain)
        target, target_mapping = read_design_chain(path, target_chain)
        if not (
            np.linalg.norm(binder.backbone[:, 1, None] - target.backbone[None, :, 1], axis=-1) <= 8
        ).any():
            raise ValueError(f"No binder-target C-alpha contacts within 8 angstroms: {path}")
        identifier = f"{source_hash[:20]}_{index}"
        reference = {
            "binder_sequence": np.array(binder.sequence),
            "target_sequence": np.array(target.sequence),
            "binder_backbone": binder.backbone,
            "target_backbone": target.backbone,
        }
        row = {
            "candidate_id": f"{campaign_id}:{epitope_id}:{arm}:{seed}:{identifier}",
            "arm": arm,
            "campaign_id": campaign_id,
            "epitope_id": epitope_id,
            "seed": seed,
            "sample_index": index,
            "binder_sequence": binder.sequence,
            "target_sequence": target.sequence,
            "source_structure": str(path),
            "source_sha256": source_hash,
            "target_reference_sha256": target_hash,
            "target_reference": str(Path(target_reference).resolve()),
            "run_record_sha256": run_hash,
            "run_record": run,
            "binder_mapping": binder_mapping,
            "target_mapping": target_mapping,
            "binder_chain": binder_chain,
            "target_chain": target_chain,
            "structure_evaluation": "pending",
            "biological_evaluation": "not_assessed",
        }
        prepared.append((identifier, reference, row))
    if len({row["target_sequence"] for _, _, row in prepared}) != 1:
        raise ValueError("An import must use one target construct")
    output.mkdir(parents=True, exist_ok=True)
    rows = []
    for identifier, reference, row in prepared:
        path = output / f"{identifier}.npz"
        np.savez_compressed(path, **reference)
        row.update(reference_npz=str(path), reference_sha256=file_hash(path))
        rows.append(row)
    (output / "candidates.jsonl").write_text("".join(json.dumps(row) + "\n" for row in rows))
    (output / "sequences.fasta").write_text(
        "".join(f">{row['candidate_id']}\n{row['binder_sequence']}\n" for row in rows)
    )
    return rows
