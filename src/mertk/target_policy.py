"""Enforce the user's active Ig1/Ig2 scope before preparing or launching new work."""

import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def load_policy():
    return json.loads((ROOT / "configs/active_target_policy.json").read_text())


def require_ig_target(epitope, hotspots, *, policy=None):
    policy = load_policy() if policy is None else policy
    domain = policy["epitopes"].get(epitope)
    if domain not in policy["allowed_domains"]:
        raise ValueError(
            f"Target {epitope!r} is retired or unreviewed; only active Ig1/Ig2 targets are allowed"
        )
    low, high = policy["allowed_domains"][domain]
    if not hotspots or any(type(n) is not int or not low <= n <= high for n in hotspots):
        raise ValueError(f"All hotspots for {epitope} must lie within {domain} A{low}-A{high}")
    return domain


def require_generation_ready(*, policy=None):
    policy = load_policy() if policy is None else policy
    if policy["new_mertk_generation_enabled"] is not True:
        raise ValueError(policy["readiness_reason"])


def require_active_target(epitope, hotspots, record, policy):
    if epitope == "Stalk_P490S491":
        scope = policy.get("cleavage_protector", {})
        if (
            scope.get("enabled") is not True
            or record.get("workflow") != "cleavage_protector_v1"
            or record.get("aim") != "aim2_shedding_protector"
        ):
            raise ValueError("Dedicated cleavage-protector authorization required")
        if hotspots != scope.get("hotspots") or any(not 485 <= n <= 498 for n in hotspots):
            raise ValueError("Protector hotspots must match the authorized cleavage-site stalk")
        return "Stalk"
    return require_ig_target(epitope, hotspots, policy=policy)


def check_manifest(manifest_path, *, experiment=None, generation=False, policy=None):
    """Check actual campaign hotspots, including manifests prepared before retirement."""
    policy = load_policy() if policy is None else policy
    manifest_path = Path(manifest_path)
    record = json.loads(manifest_path.read_text())
    mertk = False
    if "experiments" in record:
        keys = [experiment] if experiment else list(record["experiments"])
        if not keys:
            raise ValueError("No experiments specified")
        for key in keys:
            spec = record["experiments"][key]
            if key in policy["allowed_controls"] and spec["epitope"] is None:
                continue
            child = json.loads(
                (manifest_path.parent / "campaigns" / key / "manifest.json").read_text()
            )
            if child.get("epitope") != spec["epitope"]:
                raise ValueError("Suite and campaign target identities disagree")
            require_active_target(spec["epitope"], child.get("hotspots"), record, policy)
            mertk = True
    else:
        require_ig_target(record.get("epitope"), record.get("hotspots"), policy=policy)
        mertk = True
    if generation and mertk:
        require_generation_ready(policy=policy)
        if record.get("workflow") == "rfd4_native_repair_v1":
            if record.get("generation_glycan_mode") != "fixed_experimental_glycan_coordinates":
                raise ValueError("RFD4 repair requires fixed native glycan context")
            audit = json.loads((manifest_path.parent / "input_validation.json").read_text())
            for key in keys:
                spec = record["experiments"][key]
                result = audit[key]
                actual = hashlib.sha256(
                    (manifest_path.parent / spec["annotated_input"]).read_bytes()
                ).hexdigest()
                if (
                    actual != result["input_sha256"]
                    or result.get("pipeline_passed") is not True
                    or result.get("roundtrip_conditions_exact") is not True
                    or result.get("glycan_atoms", 0) < 1
                    or result.get("fixed_atoms", 0) < 1568
                ):
                    raise ValueError("RFD4 repair input validation missing or changed")
            return record
        if record.get("workflow") not in {"ig_followup_v1", "cleavage_protector_v1"}:
            raise ValueError(
                "New generation requires the glycan-guarded workflow and GPU preflight"
            )
        for key in keys:
            child_path = manifest_path.parent / "campaigns" / key
            child = json.loads((child_path / "manifest.json").read_text())
            if child.get("glycan_generation_guard") != "asn_frame_excluded_volume_v1":
                raise ValueError("Missing during-generation glycan guard")
            if not (child_path / "inputs/glycan_constraints.json").is_file():
                raise ValueError("Missing frozen glycan constraints")
            if record["workflow"] == "cleavage_protector_v1":
                spec = json.loads((child_path / "inputs/glycan_constraints.json").read_text())
                if "membrane" not in spec or spec["native_numbers"] != list(range(386, 527)):
                    raise ValueError("Protector requires membrane-aware A386–A526 context")
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--experiment")
    parser.add_argument("--generation", action="store_true")
    args = parser.parse_args()
    try:
        check_manifest(args.manifest, experiment=args.experiment, generation=args.generation)
    except ValueError as error:
        raise SystemExit(str(error)) from error
    print("ACTIVE_TARGET_POLICY_VERIFIED")


if __name__ == "__main__":
    main()
