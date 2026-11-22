"""Prepare matched, template-free AF3 jobs without importing or executing AF3."""

import hashlib
import json
import shutil
from pathlib import Path

import numpy as np

from mertk.provenance import file_hash
from mertk.structures import AA


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
    for key in ("samples_per_seed", "num_recycles"):
        if type(protocol[key]) is not int or protocol[key] < 1:
            raise ValueError(f"{key} must be a positive integer")
    if protocol["input_version"] != 1 or protocol["flash_attention_implementation"] not in {
        "triton",
        "cudnn",
        "xla",
    }:
        raise ValueError("Unsupported AF3 input version or attention implementation")
    arms = protocol["required_arms"]
    if (
        not arms
        or len(set(arms)) != len(arms)
        or any(not isinstance(a, str) or not a for a in arms)
    ):
        raise ValueError("Provide distinct nonempty required arms")
    gates = protocol["gates"]
    values = list(gates.values())
    if not all(type(x) in {int, float} and np.isfinite(x) for x in values):
        raise ValueError("Thresholds must be finite numbers")
    if not 0 <= gates["min_pair_iptm"] <= 1 or not 0 < gates["min_reference_coverage"] <= 1:
        raise ValueError("Invalid confidence or coverage threshold")
    if min(gates["max_binder_fold_rmsd"], gates["max_binder_pose_rmsd"]) <= 0:
        raise ValueError("Invalid RMSD threshold")


def validate_cohort(candidates, required_arms):
    """Compare independent designs on the same target and epitope strata."""
    if not candidates or len({c["candidate_id"] for c in candidates}) != len(candidates):
        raise ValueError("Candidate IDs must be nonempty and unique")
    arms, contexts = {}, {}
    for candidate in candidates:
        if not candidate["candidate_id"] or not candidate["arm"]:
            raise ValueError("Candidate ID and arm must be nonempty")
        if type(candidate["seed"]) is not int or candidate["seed"] < 0:
            raise ValueError("Design seed must be a nonnegative integer")
        if type(candidate["sample_index"]) is not int or candidate["sample_index"] < 0:
            raise ValueError("Sample index must be a nonnegative integer")
        context = (candidate["campaign_id"], candidate["epitope_id"])
        signature = (candidate["target_sequence"], candidate["target_reference_sha256"])
        if context in contexts and contexts[context] != signature:
            raise ValueError("Unmatched target sequence or experimental reference within a stratum")
        contexts[context] = signature
        arms.setdefault(candidate["arm"], set()).add(context)
    if not set(required_arms) <= set(arms):
        raise ValueError("Missing required comparison arms")
    if any(strata != set(contexts) for strata in arms.values()):
        raise ValueError("Every arm must cover identical campaign/epitope strata")


def checked_reference(candidate):
    if file_hash(candidate["reference_npz"]) != candidate["reference_sha256"]:
        raise ValueError("Reference structure changed after candidate export")
    with np.load(candidate["reference_npz"], allow_pickle=False) as data:
        reference = dict(data)
    binder, target = candidate["binder_sequence"], candidate["target_sequence"]
    native = str(reference["binder_sequence"])
    if str(reference["target_sequence"]) != target or native != binder:
        raise ValueError("Candidate sequence differs from its imported complex reference")
    if not binder or not target or any(x not in AA for x in binder + target):
        raise ValueError("Only standard protein residues are supported")
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
):
    validate_protocol(protocol)
    candidates = [
        json.loads(line) for line in Path(candidates_path).read_text().splitlines() if line.strip()
    ]
    validate_cohort(candidates, protocol["required_arms"])
    if samples_per_design_seed is not None:
        if samples_per_design_seed < 1:
            raise ValueError("Sample budget must be positive")
        candidates = [c for c in candidates if c["sample_index"] < samples_per_design_seed]
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
