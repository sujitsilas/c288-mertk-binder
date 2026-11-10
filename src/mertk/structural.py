"""AF3 prediction validation, fold/pose measurements and paired structural reporting."""

import csv
import json
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB.PDBExceptions import PDBConstructionException
from Bio.SeqUtils import seq1

from mertk.af3 import checked_reference, load_jobs, validate_cohort
from mertk.comparison import family_estimate
from mertk.training import file_hash


def rigid_fit(mobile, reference):
    mobile, reference = np.asarray(mobile, dtype=float), np.asarray(reference, dtype=float)
    if mobile.shape != reference.shape or mobile.ndim != 2 or mobile.shape[1] != 3:
        raise ValueError("Alignment requires matched N-by-3 coordinates")
    if len(mobile) < 3 or not np.isfinite(mobile).all() or not np.isfinite(reference).all():
        raise ValueError("Alignment requires at least three finite matched atoms")
    x, y = mobile - mobile.mean(0), reference - reference.mean(0)
    if min(np.linalg.matrix_rank(x, tol=1e-6), np.linalg.matrix_rank(y, tol=1e-6)) < 2:
        raise ValueError("Collinear coordinates cannot establish a unique pose")
    u, _, vt = np.linalg.svd(x.T @ y)
    rotation = u @ np.diag([1, 1, np.linalg.det(u @ vt)]) @ vt
    return rotation, reference.mean(0) - mobile.mean(0) @ rotation


def rmsd(a, b):
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=-1))))


def geometry_metrics(reference, target_ca, binder_ca):
    rt, rb = reference["target_backbone"][:, 1], reference["binder_backbone"][:, 1]
    if target_ca.shape != rt.shape or binder_ca.shape != rb.shape:
        raise ValueError("Predicted chains differ from full reference sequence lengths")
    mt, mb = np.isfinite(rt).all(-1), np.isfinite(rb).all(-1)
    tr, tt = rigid_fit(target_ca[mt], rt[mt])
    br, bt = rigid_fit(binder_ca[mb], rb[mb])
    ref_contacts = np.linalg.norm(rb[mb, None] - rt[None, mt], axis=-1) <= 8.0
    pred_contacts = np.linalg.norm(binder_ca[mb, None] - target_ca[None, mt], axis=-1) <= 8.0
    return {
        "binder_fold_rmsd": rmsd(binder_ca[mb] @ br + bt, rb[mb]),
        "binder_pose_rmsd": rmsd(binder_ca[mb] @ tr + tt, rb[mb]),
        "target_alignment_rmsd": rmsd(target_ca[mt] @ tr + tt, rt[mt]),
        "binder_reference_coverage": float(mb.mean()),
        "target_reference_coverage": float(mt.mean()),
        "ca_contact_retention": float(pred_contacts[ref_contacts].mean())
        if ref_contacts.any()
        else None,
    }


def read_prediction(path, target_sequence, binder_sequence):
    structure = MMCIFParser(QUIET=True, auth_chains=False, auth_residues=False).get_structure(
        "af3", str(path)
    )
    if len(structure) != 1 or {chain.id for chain in structure[0]} != {"A", "B"}:
        raise ValueError("Expected one model with explicit target A and binder B chains")
    arrays = []
    for chain_id, sequence in [("A", target_sequence), ("B", binder_sequence)]:
        ca = np.full((len(sequence), 3), np.nan)
        seen = set()
        for residue in structure[0][chain_id]:
            index = residue.id[1] - 1
            if index in seen or not 0 <= index < len(sequence) or "CA" not in residue:
                raise ValueError("Invalid, duplicate or missing predicted residue mapping")
            if seq1(residue.resname, undef_code="X") != sequence[index]:
                raise ValueError("Predicted sequence does not match the AF3 candidate")
            seen.add(index)
            ca[index] = residue["CA"].coord
        if not np.isfinite(ca).all():
            raise ValueError("Predicted structure has missing/non-finite C-alpha coordinates")
        arrays.append(ca)
    return arrays


def validate_processed_input(path, expected):
    actual = json.loads(Path(path).read_text())
    if actual["modelSeeds"] != expected["modelSeeds"]:
        raise ValueError("AF3 prediction seeds differ from the planned budget")
    chains = {}
    for entry in actual["sequences"]:
        protein = entry["protein"]
        ids = protein["id"] if isinstance(protein["id"], list) else [protein["id"]]
        for chain in ids:
            if chain in chains:
                raise ValueError("Repeated AF3 input chain ID")
            chains[chain] = protein
    if set(chains) != {"A", "B"}:
        raise ValueError("AF3 output has different input chains")
    for entry in expected["sequences"]:
        protein = entry["protein"]
        for key in ("sequence", "unpairedMsa", "pairedMsa", "templates"):
            if chains[protein["id"]].get(key) != protein[key]:
                raise ValueError(f"AF3 processed input changed {key}")


def prediction_files(root, job, protocol):
    directory = Path(root) / job["name"]
    if not (directory / f"{job['name']}_ranking_scores.csv").is_file():
        raise FileNotFoundError("Missing AF3 ranking table")
    expected_dirs = {
        f"seed-{seed}_sample-{sample}"
        for seed in protocol["model_seeds"]
        for sample in range(protocol["samples_per_seed"])
    }
    found = {p.name for p in directory.glob("seed-*_sample-*") if p.is_dir()}
    if not expected_dirs <= found:
        raise FileNotFoundError("Missing planned AF3 samples")
    if found != expected_dirs:
        raise ValueError("AF3 sample budget differs from the plan")
    samples = []
    for name in sorted(expected_dirs):
        paths = {}
        for key, suffix in [
            ("summary", "_summary_confidences.json"),
            ("confidence", "_confidences.json"),
            ("model", "_model.cif"),
        ]:
            matches = [
                p
                for p in (directory / name).glob("*" + suffix)
                if key != "confidence" or not p.name.endswith("_summary_confidences.json")
            ]
            if not matches:
                raise FileNotFoundError(f"Missing AF3 {key} for {name}")
            if len(matches) != 1:
                raise ValueError("Ambiguous AF3 prediction files")
            paths[key] = matches[0]
        samples.append((name, paths))
    return samples


def score_job(job, protocol, bundle, predictions):
    candidate = job["candidate"]
    reference = checked_reference(
        dict(candidate, reference_npz=str(Path(bundle) / candidate["reference_npz"]))
    )
    expected = json.loads((Path(bundle) / job["input"]).read_text())
    validate_processed_input(Path(predictions) / job["name"] / f"{job['name']}_data.json", expected)
    samples = prediction_files(predictions, job, protocol)
    ranking_path = Path(predictions) / job["name"] / f"{job['name']}_ranking_scores.csv"
    with ranking_path.open() as handle:
        ranking_rows = list(csv.DictReader(handle))
    ranks = {}
    for row in ranking_rows:
        name = f"seed-{int(row['seed'])}_sample-{int(row['sample'])}"
        rank = float(row["ranking_score"])
        if name in ranks or not np.isfinite(rank):
            raise ValueError("Invalid or duplicate AF3 ranking table row")
        ranks[name] = rank
    if set(ranks) != {name for name, _ in samples}:
        raise ValueError("AF3 ranking table differs from the planned sample budget")
    ranked = {}
    for name, paths in samples:
        summary = json.loads(paths["summary"].read_text())
        rank = summary["ranking_score"]
        if type(rank) not in {float, int} or not np.isfinite(rank):
            raise ValueError("Invalid AF3 ranking score")
        if abs(rank - ranks[name]) > 0.005001:
            raise ValueError("AF3 ranking table disagrees with sample confidence")
        ranked[name] = (paths, summary)
    # Summary JSON rounds scores; select using the full-precision CSV and upstream tie order.
    order = [
        f"seed-{seed}_sample-{sample}"
        for seed in protocol["model_seeds"]
        for sample in range(protocol["samples_per_seed"])
    ]
    selected = max(order, key=lambda name: ranks[name])
    paths, summary = ranked[selected]
    confidence = json.loads(paths["confidence"].read_text())
    target_ca, binder_ca = read_prediction(
        paths["model"], candidate["target_sequence"], candidate["binder_sequence"]
    )
    metrics = geometry_metrics(reference, target_ca, binder_ca)
    token_chains = np.array(confidence["token_chain_ids"])
    # Standard amino acids have one AF3 token per residue; verify full lengths.
    if (
        set(token_chains.tolist()) != {"A", "B"}
        or sum(token_chains == "A") != len(target_ca)
        or sum(token_chains == "B") != len(binder_ca)
    ):
        raise ValueError("AF3 token mapping differs from candidate chains")
    residue_ids = np.asarray(confidence["token_res_ids"])
    if residue_ids.shape != token_chains.shape:
        raise ValueError("AF3 token residue mapping has a different length")
    for chain, length in [("A", len(target_ca)), ("B", len(binder_ca))]:
        if not np.array_equal(residue_ids[token_chains == chain], np.arange(1, length + 1)):
            raise ValueError("AF3 token residue IDs differ from the full sequence mapping")
    chain_ids = summary.get("chain_ids", list(dict.fromkeys(token_chains.tolist())))
    if len(chain_ids) != 2 or set(chain_ids) != {"A", "B"}:
        raise ValueError("Invalid confidence chain mapping")
    pair = np.asarray(summary["chain_pair_iptm"], dtype=float)
    if pair.shape != (2, 2) or not np.isfinite(pair).all() or np.any((pair < 0) | (pair > 1)):
        raise ValueError("Invalid pair ipTM")
    iptm = float(pair[chain_ids.index("A"), chain_ids.index("B")])
    pae = np.asarray(confidence["pae"], dtype=float)
    if (
        pae.shape != (len(token_chains), len(token_chains))
        or not np.isfinite(pae).all()
        or (pae < 0).any()
    ):
        raise ValueError("Invalid AF3 PAE array")
    atom_chains, plddt = (
        np.array(confidence["atom_chain_ids"]),
        np.array(confidence["atom_plddts"], dtype=float),
    )
    if (
        atom_chains.shape != plddt.shape
        or plddt.ndim != 1
        or not np.isfinite(plddt).all()
        or np.any((plddt < 0) | (plddt > 100))
        or not (atom_chains == "B").any()
    ):
        raise ValueError("Invalid AF3 per-atom confidence")
    # The pinned implementation serializes this documented boolean as 0.0/1.0.
    clash = summary["has_clash"]
    if type(clash) not in {bool, int, float} or clash not in (0, 1):
        raise ValueError("AF3 clash flag must be boolean or numeric 0/1")
    clash = bool(clash)
    metrics.update(
        pair_iptm=iptm,
        has_clash=clash,
        binder_atom_plddt=float(plddt[atom_chains == "B"].mean()),
        target_to_binder_pae=float(pae[np.ix_(token_chains == "A", token_chains == "B")].mean()),
        binder_to_target_pae=float(pae[np.ix_(token_chains == "B", token_chains == "A")].mean()),
    )
    gates = protocol["gates"]
    checks = {
        "pair_confidence": iptm >= gates["min_pair_iptm"],
        "binder_fold": metrics["binder_fold_rmsd"] <= gates["max_binder_fold_rmsd"],
        "binder_pose": metrics["binder_pose_rmsd"] <= gates["max_binder_pose_rmsd"],
        "reference_coverage": min(
            metrics["target_reference_coverage"], metrics["binder_reference_coverage"]
        )
        >= gates["min_reference_coverage"],
        "no_severe_clash": not clash,
    }
    return {
        **candidate,
        "prediction_status": "complete",
        "structure_evaluation": "scored",
        "structural_pass": all(checks.values()),
        "gates": checks,
        "metrics": metrics,
        "selected_sample": selected,
        "predictions": len(samples),
        "ranking_score": ranks[selected],
        "ranking_table_sha256": file_hash(ranking_path),
        "selected_files": {
            key: {"path": str(path), "sha256": file_hash(path)} for key, path in paths.items()
        },
        "biological_evaluation": "not_assessed",
    }


def summarize_structures(records, protocol, execution_verified=False):
    validate_cohort(records, protocol["required_arms"])
    complete = all(row["prediction_status"] == "complete" for row in records)
    arms = {}
    for row in records:
        arms.setdefault(row["arm"], []).append(row)
    report = {
        "endpoint": "predicted_structural_consistency",
        "protocol_state": protocol["protocol_state"],
        "complete_cohort": complete,
        "execution_verified": execution_verified,
        "biological_evaluation": "not_assessed",
        "arms": {},
    }
    baselines = {
        name: {(row["complex_id"], row["seed"], row["sample_index"]): row for row in arms[name]}
        for name in ("proteinmpnn", "lora", "frozen_esm3")
        if name in arms
    }
    for arm, rows in arms.items():
        valid = [row for row in rows if row["prediction_status"] == "complete"]
        estimate = {
            "planned_candidates": len(rows),
            "valid_predictions": len(valid),
            "valid_output_rate": len(valid) / len(rows),
            "passing_candidates": sum(row["structural_pass"] for row in valid),
        }
        # No complete-case deletion: incomplete cohorts do not receive comparative estimates.
        if complete:
            grouped = {}
            for row in rows:
                grouped.setdefault(row["component"], []).append(row)
            rates = [
                np.mean([row["structural_pass"] for row in group]) for group in grouped.values()
            ]
            estimate["family_pass_rate"] = family_estimate(rates, np.random.default_rng(42))
            if not rows[0]["control"]:
                for baseline_name, baseline in baselines.items():
                    differences = [
                        np.mean(
                            [
                                int(row["structural_pass"])
                                - int(
                                    baseline[(row["complex_id"], row["seed"], row["sample_index"])][
                                        "structural_pass"
                                    ]
                                )
                                for row in group
                            ]
                        )
                        for group in grouped.values()
                    ]
                    estimate[f"paired_pass_rate_delta_vs_{baseline_name}"] = family_estimate(
                        differences, np.random.default_rng(42)
                    )
                paired = estimate["paired_pass_rate_delta_vs_proteinmpnn"]
                eligible = (
                    arm != "proteinmpnn"
                    and execution_verified
                    and protocol["protocol_state"] == "locked"
                    and rows[0]["split"] == "test"
                    and len(grouped) >= max(2, protocol["min_families_for_inference"])
                )
                estimate["noninferiority_supported"] = (
                    bool(paired["ci_95"][0] > -protocol["noninferiority_margin"])
                    if eligible
                    else None
                )
        report["arms"][arm] = estimate
    return report


def score_bundle(manifest_path, predictions, output):
    manifest_path, predictions, output = (
        Path(manifest_path).resolve(),
        Path(predictions).resolve(),
        Path(output),
    )
    manifest = load_jobs(manifest_path)
    if output.exists() and any(output.iterdir()):
        raise FileExistsError("Choose an empty structural report directory")
    records = []
    for job in manifest["jobs"]:
        try:
            record = score_job(job, manifest["protocol"], manifest_path.parent, predictions)
        except (
            FileNotFoundError,
            ValueError,
            KeyError,
            TypeError,
            IndexError,
            PDBConstructionException,
        ) as error:
            record = {
                **job["candidate"],
                "prediction_status": "missing_output"
                if isinstance(error, FileNotFoundError)
                else "invalid_output",
                "error": str(error),
                "structural_pass": None,
                "structure_evaluation": "incomplete",
                "biological_evaluation": "not_assessed",
            }
        records.append(record)
    execution_file = predictions / "execution.json"
    execution = json.loads(execution_file.read_text()) if execution_file.exists() else None
    verified = bool(
        execution
        and execution["plan_id"] == manifest["plan_id"]
        and execution["engine"]["revision"] == manifest["engine"]["revision"]
        and execution["engine"].get("weights")
        and execution["engine"].get("environment")
        and {job["name"] for job in manifest["jobs"]}
        <= {name for attempt in execution.get("attempts", []) for name in attempt["candidates"]}
    )
    report = summarize_structures(records, manifest["protocol"], verified)
    report.update(plan_id=manifest["plan_id"], protocol=manifest["protocol"], execution=execution)
    output.mkdir(parents=True, exist_ok=True)
    (output / "candidates.jsonl").write_text("".join(json.dumps(row) + "\n" for row in records))
    (output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    return report
