"""Matched sequence candidates and immutable references for structural evaluation."""

import json
from pathlib import Path

import numpy as np

from mertk.comparison import sequence_metrics, summarize
from mertk.training import file_hash


def export_candidates(rows, results, output):
    summarize(results)  # Fail before exporting unmatched arm/seed cohorts.
    if any(result["smoke_test"] for result in results):
        raise ValueError("Synthetic training results cannot become structural benchmark candidates")
    application = results[0].get("reference_kind") == "bindcraft_seed"
    row_index = {row["id"]: row for row in rows}
    contexts, examples, candidates = {}, {}, []
    for case in results[0]["cases"]:
        row = row_index[case["id"]]
        if row["component"] != case["component"] or row["split"] != results[0]["split"]:
            raise ValueError("Candidate family or split changed after comparison")
        if file_hash(row["path"]) != case["source_sha256"]:
            raise ValueError("Prepared structure changed after sequence comparison")
        if "source_pdb" in row and file_hash(row["source_pdb"]) != row["source_sha256"]:
            raise ValueError("Original complex changed after import")
        with np.load(row["path"], allow_pickle=False) as data:
            example = dict(data)
        examples[row["id"]] = example
        context = {
            "complex_id": row["id"],
            "component": row["component"],
            "split": row["split"],
            "reference_npz": str(Path(row["path"]).resolve()),
            "reference_sha256": case["source_sha256"],
            "mask_sha256": case["mask_sha256"],
            "dataset_identity": results[0]["dataset_identity"],
            "target_sequence": str(example["target_sequence"]),
            "design_positions_1based": (np.flatnonzero(example["design_mask"]) + 1).tolist(),
            "reference_kind": results[0]["reference_kind"],
            "smoke_test": False,
            "structure_evaluation": "pending",
            "biological_evaluation": "not_assessed",
        }
        contexts[row["id"]] = context
        control = "bindcraft_original" if application else "native_original"
        candidates.append(
            {
                **context,
                "candidate_id": f"{row['id']}:{control}",
                "arm": control,
                "control": True,
                "seed": None,
                "sample_index": 0,
                "binder_sequence": str(example["binder_sequence"]),
            }
        )
    for result in results:
        for case in result["cases"]:
            example = examples[case["id"]]
            sequence_metrics(
                str(example["binder_sequence"]), case["sequences"], example["design_mask"]
            )
            for index, sequence in enumerate(case["sequences"]):
                candidates.append(
                    {
                        **contexts[case["id"]],
                        "candidate_id": f"{case['id']}:{result['arm']}:{result['seed']}:{index}",
                        "arm": result["arm"],
                        "control": False,
                        "seed": result["seed"],
                        "sample_index": index,
                        "binder_sequence": sequence,
                        "checkpoint": result.get("checkpoint"),
                        "checkpoint_sha256": result.get("checkpoint_sha256"),
                        "model_identity": result["model_identity"],
                    }
                )
    Path(output).write_text("".join(json.dumps(candidate) + "\n" for candidate in candidates))
    return candidates
