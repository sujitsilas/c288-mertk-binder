"""Matched sequence-design metrics and family-level paired comparisons."""

import hashlib
import json
import time
from itertools import combinations
from pathlib import Path

import numpy as np
import torch

from mertk.data import AA
from mertk.dataset import dataset_fingerprint, make_batch, read_records, synthetic_example
from mertk.model import weight_identity
from mertk.objectives import sample_sequence
from mertk.runtime import autocast, resolve_device, resolve_precision, synchronize
from mertk.training import audit_splits, build_model, file_hash


def sequence_metrics(reference, sequences, mask):
    mask = np.asarray(mask, dtype=bool)
    if not mask.any() or len(mask) != len(reference) or len(sequences) < 2:
        raise ValueError("Metrics require a nonempty aligned mask and at least two samples")
    reference = np.array(list(reference))
    arrays = []
    for sequence in sequences:
        array = np.array(list(sequence))
        if len(array) != len(reference) or not np.array_equal(array[~mask], reference[~mask]):
            raise ValueError("Generated sequence changed fixed positions or length")
        if any(aa not in AA for aa in array[mask]):
            raise ValueError("Designed positions must be standard amino acids")
        arrays.append(array[mask])
    return {
        "mean_recovery": float(np.mean([np.mean(a == reference[mask]) for a in arrays])),
        "pairwise_diversity": float(np.mean([np.mean(a != b) for a, b in combinations(arrays, 2)])),
        "unique_fraction": len(set(sequences)) / len(sequences),
        "design_residues": int(mask.sum()),
        "samples": len(sequences),
    }


def dataset_identity(manifest, smoke=False):
    if smoke:
        return "synthetic-smoke-v2"
    return dataset_fingerprint(manifest)


def run_comparison(
    config,
    arm,
    output,
    checkpoint=None,
    samples=8,
    sampling_steps=16,
    temperature=1.0,
    split="test",
    smoke=False,
    expected_dataset=None,
    application_manifest=None,
):
    output = Path(output)
    if output.exists():
        raise FileExistsError(output)
    if samples < 2:
        raise ValueError("Use at least two samples to measure diversity")
    manifest = application_manifest or config["manifest"]
    if application_manifest and (split != "application_holdout" or smoke):
        raise ValueError("Application comparison requires real application_holdout examples")
    identity = dataset_identity(manifest, smoke)
    if expected_dataset is not None and expected_dataset != identity:
        raise ValueError("Dataset changed since the suite was planned")
    device = resolve_device("cpu" if smoke else config["device"])
    torch.manual_seed(config["seed"])
    if config.get("cpu_threads"):
        torch.set_num_threads(config["cpu_threads"])
    started = time.perf_counter()
    training_subset = None
    if smoke:
        rows = [{"id": "synthetic", "component": "synthetic", "split": "test"}]
    else:
        rows = read_records(manifest)
        audit_splits(rows)
        if application_manifest:
            if any(row["split"] != "application_holdout" for row in rows):
                raise ValueError("Application manifest must contain only holdouts")
            audit_splits(read_records(config["manifest"]) + rows)
        rows = [row for row in rows if row["split"] == split]
        if not rows:
            raise ValueError(f"No examples for {split}")
    if arm == "proteinmpnn":
        from mertk.proteinmpnn import load_proteinmpnn

        model, model_identity = load_proteinmpnn(device)
    else:
        model = build_model(config, device, tiny=smoke).eval()
        model_identity = (
            {"kind": "random-tiny-smoke"}
            if smoke
            else (
                {"kind": "atom_decoder"}
                if config.get("model_type") == "atom_decoder"
                else weight_identity(config["weights"])
            )
        )
        if checkpoint:
            saved = torch.load(checkpoint, map_location="cpu", weights_only=True)
            if (
                saved["metadata"]["smoke_test"] != smoke
                or saved["metadata"]["weights"] != model_identity
            ):
                raise ValueError("Checkpoint type or base weights mismatch")
            for key in (
                "model_type",
                "atom_context",
                "objective",
                "use_lora",
                "use_geometry",
                "rank",
                "alpha",
                "hidden",
                "seed",
                "sampling",
                "data_seed",
                "train_family_fraction",
            ):
                if saved["metadata"]["config"].get(key) != config.get(key):
                    raise ValueError(f"Checkpoint arm mismatch: {key}")
            if not smoke and saved["metadata"]["manifest_sha256"] != file_hash(config["manifest"]):
                raise ValueError("Checkpoint was trained on a different split manifest")
            if (
                not smoke
                and "dataset_identity" in saved["metadata"]
                and saved["metadata"]["dataset_identity"] != dataset_identity(config["manifest"])
            ):
                raise ValueError("Prepared training data differ from checkpoint inputs")
            model.load_adapter_state_dict(saved["adapter"])
            if "training_examples" in saved["metadata"]:
                training_subset = hashlib.sha256(
                    json.dumps(saved["metadata"]["training_examples"]).encode()
                ).hexdigest()
        elif arm != "frozen_esm3":
            raise ValueError("Trainable comparison arms require a checkpoint")
    synchronize(device)
    load_seconds = time.perf_counter() - started
    cases = []
    for row in rows:
        if smoke:
            example = synthetic_example()
            source_hash = identity
        else:
            with np.load(row["path"], allow_pickle=False) as data:
                example = dict(data)
            source_hash = file_hash(row["path"])
        reference = str(example["binder_sequence"])
        mask = example["design_mask"]
        sequences = []
        if arm != "proteinmpnn":
            batch = make_batch(example, model.tokenizer, device, require_structure=not smoke)
        synchronize(device)
        started = time.perf_counter()
        for i in range(samples):
            # Stable per-example seed: unchanged by row order, model size, or prior examples.
            seed = int(
                hashlib.sha256(f"{config['seed']}:{row['id']}:{i}".encode()).hexdigest()[:8], 16
            )
            if arm == "proteinmpnn":
                from mertk.proteinmpnn import sample_proteinmpnn

                sequence = sample_proteinmpnn(model, example, device, seed, temperature)
            else:
                with autocast(device, "fp32" if smoke else config.get("precision", "auto")):
                    tokens = sample_sequence(
                        model, batch, model.tokenizer, sampling_steps, temperature, seed
                    )
                fixed = ~batch["design_mask"]
                if not torch.equal(tokens[fixed], batch["sequence_tokens"][fixed]):
                    raise RuntimeError("Sampler changed fixed target/context tokens")
                sequence = "".join(
                    model.tokenizer.convert_ids_to_tokens(
                        tokens[0, 1 : len(reference) + 1].tolist()
                    )
                )
            sequences.append(sequence)
        synchronize(device)
        sampling_seconds = time.perf_counter() - started
        metrics = sequence_metrics(reference, sequences, mask)
        if application_manifest:
            metrics["mean_seed_retention"] = metrics.pop("mean_recovery")
        cases.append(
            {
                "id": row["id"],
                "component": row["component"],
                "source_sha256": source_hash,
                "mask_sha256": hashlib.sha256(mask.astype(np.uint8).tobytes()).hexdigest(),
                "seed_sequence": reference,
                "sequences": sequences,
                "sampling_seconds": sampling_seconds,
                **metrics,
            }
        )
    result = {
        "arm": arm,
        "seed": config["seed"],
        "split": split,
        "smoke_test": smoke,
        "reference_kind": "bindcraft_seed" if application_manifest else "native",
        "dataset_identity": identity,
        "model_identity": model_identity,
        "sampling": {"samples": samples, "steps": sampling_steps, "temperature": temperature},
        "sampling_protocol": config.get("sampling_protocol"),
        "checkpoint": str(checkpoint) if checkpoint else None,
        "checkpoint_sha256": file_hash(checkpoint) if checkpoint else None,
        "training_subset_sha256": training_subset,
        "runtime": {
            "device": str(device),
            "torch": str(torch.__version__),
            "precision": "fp32"
            if smoke or arm == "proteinmpnn"
            else resolve_precision(device, config.get("precision", "auto")),
            "model_load_seconds": load_seconds,
            "sampling_seconds": sum(case["sampling_seconds"] for case in cases),
        },
        "cases": cases,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n")
    return result


def family_estimate(values, rng):
    values = np.array(values, dtype=float)
    ci = None
    if len(values) >= 2:
        means = np.mean(rng.choice(values, size=(2000, len(values)), replace=True), axis=1)
        ci = np.quantile(means, [0.025, 0.975]).tolist()
    return {"mean": float(values.mean()), "ci_95": ci, "families": len(values)}


def summarize(results):
    if not results:
        raise ValueError("No comparison results")
    for result in results:
        protocol = result.get("sampling_protocol")
        if protocol and result["sampling"] != {
            "samples": protocol["samples"],
            "steps": protocol["steps"],
            "temperature": protocol["temperatures"][result["arm"]],
        }:
            raise ValueError("Result differs from its declared sampling protocol")
    context = {
        (
            r["dataset_identity"],
            r["split"],
            r["smoke_test"],
            r.get("reference_kind", "native"),
            json.dumps(r.get("sampling_protocol") or r["sampling"], sort_keys=True),
        )
        for r in results
    }
    if len(context) != 1:
        raise ValueError("Comparison results have different datasets, splits, or sampling settings")
    subsets = {r["training_subset_sha256"] for r in results if r.get("training_subset_sha256")}
    if len(subsets) > 1:
        raise ValueError("Trainable arms used different training family subsets")
    arms = {}
    for result in results:
        arm = arms.setdefault(result["arm"], {})
        for case in result["cases"]:
            key = (case["id"], result["seed"])
            if key in arm:
                raise ValueError("Duplicate example/seed in comparison")
            arm[key] = case
    if "frozen_esm3" not in arms:
        raise ValueError("The matched frozen ESM3 baseline is required")
    reference = arms["frozen_esm3"]
    report = {
        "dataset_identity": results[0]["dataset_identity"],
        "split": results[0]["split"],
        "smoke_test": results[0]["smoke_test"],
        "reference_kind": results[0].get("reference_kind", "native"),
        "arms": {},
    }
    reference_metric = (
        "mean_seed_retention" if report["reference_kind"] == "bindcraft_seed" else "mean_recovery"
    )
    delta_prefix = (
        "paired_seed_retention_delta_vs_"
        if report["reference_kind"] == "bindcraft_seed"
        else "paired_recovery_delta_vs_"
    )
    # Validate every arm before computing any paired contrasts, including MPNN.
    for name, cases in arms.items():
        if set(cases) != set(reference):
            raise ValueError(f"Missing matched examples/seeds in {name}")
        for key, case in cases.items():
            base = reference[key]
            if any(
                case[field] != base[field]
                for field in ("source_sha256", "mask_sha256", "component", "design_residues")
            ):
                raise ValueError(f"Unmatched inputs or masks in {name}: {key}")
    baselines = {"frozen": reference}
    if "proteinmpnn" in arms:
        baselines["proteinmpnn"] = arms["proteinmpnn"]
    if "lora" in arms:
        baselines["lora"] = arms["lora"]
    for name, cases in arms.items():
        grouped = {}
        for key, case in cases.items():
            grouped.setdefault(case["component"], []).append(key)
        estimates = {}
        for metric in (reference_metric, "pairwise_diversity", "unique_fraction"):
            values = [np.mean([cases[key][metric] for key in keys]) for keys in grouped.values()]
            estimates[metric] = family_estimate(values, np.random.default_rng(42))
        for baseline_name, baseline in baselines.items():
            differences = [
                np.mean(
                    [cases[key][reference_metric] - baseline[key][reference_metric] for key in keys]
                )
                for keys in grouped.values()
            ]
            estimates[delta_prefix + baseline_name] = family_estimate(
                differences, np.random.default_rng(42)
            )
        report["arms"][name] = {"example_seed_pairs": len(cases), **estimates}
    return report
