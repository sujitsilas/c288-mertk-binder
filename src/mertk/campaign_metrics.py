"""Count independent trajectories separately from variants and diagnostic files."""

import csv
import json
import math
import re
from pathlib import Path

from mertk.provenance import file_hash


def csv_rows(path):
    if not Path(path).is_file():
        return []
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle))


def acceptance_interval(successes, attempts):
    """Wilson two-sided 95% CI; also exact one-sided 95% bound for zero successes."""
    if attempts < 0 or not 0 <= successes <= attempts:
        raise ValueError("Successes must lie between zero and attempts")
    if attempts == 0:
        return {"rate": None, "wilson_95": None, "zero_hit_upper_95_one_sided": None}
    z = 1.959963984540054
    rate = successes / attempts
    denom = 1 + z * z / attempts
    centre = (rate + z * z / (2 * attempts)) / denom
    width = z * math.sqrt(rate * (1 - rate) / attempts + z * z / (4 * attempts**2)) / denom
    return {
        "rate": rate,
        "wilson_95": [max(0, centre - width), min(1, centre + width)],
        "zero_hit_upper_95_one_sided": 1 - 0.05 ** (1 / attempts) if successes == 0 else None,
    }


def summarize(campaign, *, log_path=None, completed=None):
    campaign = Path(campaign)
    manifest = json.loads((campaign / "manifest.json").read_text())
    designs = campaign / "bindcraft/designs"
    log = Path(log_path) if log_path is not None else campaign / "bindcraft/run.log"
    if log_path is None and not log.is_file():
        # The native log is the direct capture. Scheduler stdout can duplicate that log.
        scheduler_logs = sorted(
            (campaign / "logs").glob(f"{manifest['campaign_id']}_bindcraft_generate_*.log")
        )
        if len(scheduler_logs) > 1:
            raise ValueError(
                "Multiple generation logs; select one explicitly to count trajectories"
            )
        if scheduler_logs:
            log = scheduler_logs[0]
    if not log.is_file():
        raise FileNotFoundError(f"No generation log for {manifest['campaign_id']}: {log}")
    text = log.read_text(errors="replace")
    if not text.strip():
        raise ValueError(f"Generation log is empty: {log}")
    starts = re.findall(r"^Starting trajectory: (\S+)", text, re.M)
    if len(starts) != len(set(starts)):
        raise ValueError("Duplicate trajectory starts in one run")
    if completed is None:
        completed = "PROJECT_RUN_COMPLETED" in text or "Finished all designs." in text
    trajectories = csv_rows(designs / "trajectory_stats.csv")
    final = csv_rows(designs / "final_design_stats.csv")
    accepted_files = sorted((designs / "Accepted").glob("*.pdb"))
    accepted_names = {re.sub(r"_model\d+$", "", p.stem) for p in accepted_files}
    final_names = {r["Design"] for r in final}
    accepted = accepted_names & final_names
    parents = {re.sub(r"_mpnn\d+$", "", name) for name in accepted}
    if not parents <= set(starts):
        raise ValueError("Accepted trajectory missing from the generation log")
    # If interrupted, the last admitted trajectory may have incomplete MPNN validation.
    finished = starts if completed else starts[:-1]
    finished_successes = len(parents & set(finished))
    failures = csv_rows(designs / "failure_csv.csv")
    sources = [p for p in designs.glob("*.csv") if p.is_file()] + [log]
    quality_path = campaign / "native_quality.json"
    quality = json.loads(quality_path.read_text()) if quality_path.exists() else {}
    return {
        "campaign_id": manifest["campaign_id"],
        "epitope": manifest.get("epitope"),
        "attempts_started": len(starts),
        "attempts_completed": len(finished),
        "attempts_censored": len(starts) - len(finished),
        "run_completed_normally": completed,
        "trajectories_scored": len(trajectories),
        "saved_sequence_proposals": len(list((designs / "MPNN/Sequences").glob("*.fasta"))),
        "sequences_with_full_metrics": len(csv_rows(designs / "mpnn_design_stats.csv")),
        "base_af2_rejected_variants": len(re.findall(r"Base AF2 filters not passed for", text)),
        "accepted_variants": len(accepted),
        "accepted_trajectories": len(parents),
        "accepted_evidence_mismatch": sorted(accepted_names ^ final_names),
        "acceptance_on_completed_attempts": acceptance_interval(finished_successes, len(finished)),
        "trajectory_names": starts,
        "accepted_names": sorted(accepted),
        "failure_counts_overlapping": {
            k: int(float(v))
            for k, v in (failures[0] if failures else {}).items()
            if v and float(v) != 0
        },
        "trajectory_metrics": [
            {
                k: row[k]
                for k in [
                    "Design",
                    "pLDDT",
                    "i_pTM",
                    "i_pAE",
                    "Surface_Hydrophobicity",
                    "Target_RMSD",
                ]
                if k in row
            }
            for row in trajectories
        ],
        "native_context": {
            arm: {k: v for k, v in data.items() if k != "records"}
            for arm, data in quality.get("arms", {}).items()
        },
        "source_sha256": {str(p): file_hash(p) for p in sources},
        "limitations": [
            "Saved FASTAs count proposed variants, not independent trajectories or all raw MPNN draws.",
            "Failure counters overlap. Trajectory pLDDT/iPAE are upstream normalized values, not AF3 units.",
            "Intervals are descriptive binomial summaries of completed attempts under these settings; "
            "do not pool different epitopes or compare biased initial-guess and default runs as equivalent.",
            "Unfinished last attempts are censored; a running or time-limited run is not a completed fixed-N experiment.",
        ],
    }
