"""Prepare matched, template-free AF3 jobs without importing or executing AF3."""

import hashlib
import json
import shutil
from pathlib import Path

import numpy as np

from mertk.data import AA
from mertk.training import file_hash


def object_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def validate_protocol(protocol):
    if protocol["protocol_state"] not in {"draft", "locked"}:
        raise ValueError("Protocol state must be draft or locked")
    if protocol["target_msa"] not in {"provided", "single_sequence"}:
        raise ValueError("Target MSA must be provided or explicitly single_sequence")
    seeds = protocol["model_seeds"]
    if (
        not seeds
        or len(set(seeds)) != len(seeds)
        or any(type(s) is not int or s < 0 for s in seeds)
    ):
        raise ValueError("AF3 seeds must be distinct nonnegative integers")
    for key in ("samples_per_seed", "num_recycles", "min_families_for_inference"):
        if type(protocol[key]) is not int or protocol[key] < 1:
            raise ValueError(f"{key} must be a positive integer")
    if protocol["input_version"] != 1 or protocol["flash_attention_implementation"] not in {
        "triton",
        "cudnn",
        "xla",
    }:
        raise ValueError("Unsupported AF3 input version or attention implementation")
    if "proteinmpnn" not in protocol["required_arms"] or len(set(protocol["required_arms"])) != len(
        protocol["required_arms"]
    ):
        raise ValueError("Distinct required arms must include proteinmpnn")
    gates = protocol["gates"]
    values = [*gates.values(), protocol["noninferiority_margin"]]
    if not all(type(x) in {int, float} and np.isfinite(x) for x in values):
        raise ValueError("Thresholds must be finite numbers")
    if not 0 <= gates["min_pair_iptm"] <= 1 or not 0 < gates["min_reference_coverage"] <= 1:
        raise ValueError("Invalid confidence or coverage threshold")
    if (
        min(gates["max_binder_fold_rmsd"], gates["max_binder_pose_rmsd"]) <= 0
        or not 0 < protocol["noninferiority_margin"] < 1
    ):
        raise ValueError("Invalid RMSD threshold or noninferiority margin")


def validate_cohort(candidates, required_arms):
    if not candidates or len({c["candidate_id"] for c in candidates}) != len(candidates):
        raise ValueError("Candidate IDs must be nonempty and unique")
    arms, references, controls = {}, {}, set()
    contexts = set()
    for candidate in candidates:
        if candidate.get("smoke_test", True):
            raise ValueError("AF3 benchmarking requires real, non-smoke candidates")
        if candidate["split"] not in {"val", "test", "application_holdout"}:
            raise ValueError("Cannot evaluate training examples as a structural benchmark")
        if (
            candidate["reference_kind"] not in {"native", "bindcraft_seed"}
            or type(candidate["control"]) is not bool
        ):
            raise ValueError("Invalid candidate reference or control type")
        if candidate["control"]:
            expected_control = (
                "native_original"
                if candidate["reference_kind"] == "native"
                else "bindcraft_original"
            )
            if (
                candidate["arm"] != expected_control
                or candidate["seed"] is not None
                or candidate["sample_index"] != 0
            ):
                raise ValueError("Original controls need their own arm and no design seed")
        elif any(
            type(candidate[key]) is not int or candidate[key] < 0
            for key in ("seed", "sample_index")
        ):
            raise ValueError("Design seed and sample index must be nonnegative integers")
        contexts.add(
            (candidate["dataset_identity"], candidate["split"], candidate["reference_kind"])
        )
        signature = tuple(
            candidate[key]
            for key in ("component", "reference_sha256", "mask_sha256", "target_sequence")
        )
        complex_id = candidate["complex_id"]
        if complex_id in references and references[complex_id] != signature:
            raise ValueError("Unmatched target, structure, family or design mask")
        references[complex_id] = signature
        if candidate["control"]:
            if complex_id in controls:
                raise ValueError("Only one original control is allowed per complex")
            controls.add(complex_id)
            continue
        key = (complex_id, candidate["seed"], candidate["sample_index"])
        arm = arms.setdefault(candidate["arm"], set())
        if key in arm:
            raise ValueError("Duplicate complex/seed/sample within an arm")
        arm.add(key)
    if len(contexts) != 1 or not set(required_arms) <= set(arms):
        raise ValueError("Mixed datasets or missing required comparison arms")
    reference = arms["proteinmpnn"]
    if controls != {key[0] for key in reference}:
        raise ValueError("Each complex needs its original sequence control")
    if any(keys != reference for keys in arms.values()):
        raise ValueError("Every arm must have identical complex/seed/sample coverage")


def checked_reference(candidate):
    if file_hash(candidate["reference_npz"]) != candidate["reference_sha256"]:
        raise ValueError("Reference structure changed after candidate export")
    with np.load(candidate["reference_npz"], allow_pickle=False) as data:
        reference = dict(data)
    binder, target = candidate["binder_sequence"], candidate["target_sequence"]
    native = str(reference["binder_sequence"])
    mask = reference["design_mask"].astype(bool)
    if hashlib.sha256(mask.astype(np.uint8).tobytes()).hexdigest() != candidate["mask_sha256"]:
        raise ValueError("Reference design mask changed")
    if str(reference["target_sequence"]) != target or len(native) != len(binder):
        raise ValueError("Candidate sequence length or target differs from reference")
    if not binder or not target or any(x not in AA + "X" for x in binder + target):
        raise ValueError("Only standard protein residues and X are supported")
    if candidate["control"]:
        if binder != native:
            raise ValueError("An original control must retain its original sequence")
    elif any(a != b and not design for a, b, design in zip(binder, native, mask, strict=True)):
        raise ValueError("Candidate changed a fixed binder residue")
    if not candidate["control"] and any(
        aa not in AA for aa, design in zip(binder, mask, strict=True) if design
    ):
        raise ValueError("Designed positions must be standard amino acids")
    if (np.flatnonzero(mask) + 1).tolist() != candidate["design_positions_1based"]:
        raise ValueError("Candidate design positions differ from reference")
    return reference


def read_target_msa(path, sequence):
    text = Path(path).read_text()
    records = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            records.append("")
        elif records:
            records[-1] += line
        else:
            raise ValueError("Target MSA must be A3M/FASTA")
    if not records or records[0] != sequence:
        raise ValueError("Target MSA query must exactly match the target sequence")
    for record in records:
        aligned = "".join(x for x in record if not x.islower())
        if len(aligned) != len(sequence) or any(x not in AA + "X-" for x in aligned):
            raise ValueError("Invalid aligned target MSA row")
    return text


def prepare_jobs(
    candidates_path,
    protocol,
    output,
    target_msas=None,
    samples_per_design_seed=None,
    max_complexes=None,
):
    validate_protocol(protocol)
    candidates = [
        json.loads(line) for line in Path(candidates_path).read_text().splitlines() if line.strip()
    ]
    validate_cohort(candidates, protocol["required_arms"])
    if samples_per_design_seed is not None:
        if samples_per_design_seed < 1:
            raise ValueError("Sample budget must be positive")
        candidates = [
            c for c in candidates if c["control"] or c["sample_index"] < samples_per_design_seed
        ]
    if max_complexes is not None:
        if max_complexes < 1:
            raise ValueError("Complex budget must be positive")
        chosen = set(
            sorted({c["complex_id"] for c in candidates}, key=lambda x: object_hash(x))[
                :max_complexes
            ]
        )
        candidates = [c for c in candidates if c["complex_id"] in chosen]
    validate_cohort(candidates, protocol["required_arms"])
    output = Path(output).resolve()
    if output.exists() and any(output.iterdir()):
        raise FileExistsError("Choose an empty AF3 job directory")
    msas = {}
    for candidate in candidates:
        checked_reference(candidate)
        sequence = candidate["target_sequence"]
        key = hashlib.sha256(sequence.encode()).hexdigest()
        if key not in msas:
            if protocol["target_msa"] == "provided":
                if not target_msas or key not in target_msas:
                    raise ValueError(f"Provide a cached target MSA for sequence SHA256 {key}")
                msas[key] = read_target_msa(target_msas[key], sequence)
            else:
                msas[key] = ""
    (output / "inputs").mkdir(parents=True)
    (output / "references").mkdir()
    jobs = []
    for candidate in candidates:
        name = "c" + object_hash(candidate["candidate_id"])[:24]
        target_hash = hashlib.sha256(candidate["target_sequence"].encode()).hexdigest()
        job_input = {
            "name": name,
            "modelSeeds": protocol["model_seeds"],
            "dialect": "alphafold3",
            "version": 1,
            "sequences": [
                {
                    "protein": {
                        "id": "A",
                        "sequence": candidate["target_sequence"],
                        "unpairedMsa": msas[target_hash],
                        "pairedMsa": "",
                        "templates": [],
                    }
                },
                {
                    "protein": {
                        "id": "B",
                        "sequence": candidate["binder_sequence"],
                        "unpairedMsa": "",
                        "pairedMsa": "",
                        "templates": [],
                    }
                },
            ],
        }
        input_path = Path("inputs") / f"{name}.json"
        (output / input_path).write_text(json.dumps(job_input, indent=2) + "\n")
        reference_path = Path("references") / f"{candidate['reference_sha256']}.npz"
        if not (output / reference_path).exists():
            shutil.copy2(candidate["reference_npz"], output / reference_path)
        jobs.append(
            {
                "name": name,
                "input": str(input_path),
                "input_sha256": file_hash(output / input_path),
                "candidate": dict(candidate, reference_npz=str(reference_path)),
            }
        )
    definition = {
        "protocol": protocol,
        "candidates_sha256": file_hash(candidates_path),
        "jobs": jobs,
        "engine": json.loads(Path("configs/external_tools.json").read_text())["alphafold3"],
    }
    manifest = {"plan_id": object_hash(definition), **definition}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def load_jobs(path):
    path = Path(path).resolve()
    manifest = json.loads(path.read_text())
    if manifest["plan_id"] != object_hash({k: v for k, v in manifest.items() if k != "plan_id"}):
        raise ValueError("AF3 plan changed after preparation")
    validate_protocol(manifest["protocol"])
    for job in manifest["jobs"]:
        if file_hash(path.parent / job["input"]) != job["input_sha256"]:
            raise ValueError("AF3 job input changed after preparation")
        checked_reference(
            dict(
                job["candidate"], reference_npz=str(path.parent / job["candidate"]["reference_npz"])
            )
        )
    validate_cohort(
        [job["candidate"] for job in manifest["jobs"]], manifest["protocol"]["required_arms"]
    )
    return manifest
