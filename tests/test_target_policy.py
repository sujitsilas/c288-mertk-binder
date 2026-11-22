"""Retired targets and unvalidated glycan handling must not consume new GPU time."""

import json

import pytest

from mertk.target_policy import check_manifest, require_generation_ready, require_ig_target


@pytest.mark.parametrize(
    "epitope,hotspots", [("FNIIIb_F472", [469, 471, 472]), ("patchA", [365, 367, 369]), (None, [])]
)
def test_retired_or_unrestricted_target_rejected(epitope, hotspots):
    with pytest.raises(ValueError, match="retired or unreviewed"):
        require_ig_target(epitope, hotspots)


def test_ig_name_cannot_hide_non_ig_hotspots():
    with pytest.raises(ValueError, match="All hotspots"):
        require_ig_target("Ig2_L270", [283, 308, 367])


def test_preparation_allows_both_ig_domains():
    assert require_ig_target("Ig1_E152", [122, 148]) == "Ig1"
    assert require_ig_target("Ig2_L270", [196, 270]) == "Ig2"


def test_generation_waits_for_real_glycan_implementation():
    with pytest.raises(ValueError, match="glycan exclusion"):
        require_generation_ready(
            policy={
                "new_mertk_generation_enabled": False,
                "readiness_reason": "During-generation glycan exclusion not validated",
            }
        )


def test_historical_round_cannot_be_relaunched(tmp_path):
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({"experiments": {"FNIIIb_F472": {"epitope": "FNIIIb_F472"}}}))
    child = tmp_path / "campaigns/FNIIIb_F472/manifest.json"
    child.parent.mkdir(parents=True)
    child.write_text(json.dumps({"epitope": "FNIIIb_F472", "hotspots": [469, 471, 472]}))
    with pytest.raises(ValueError, match="FNIIIb_F472"):
        check_manifest(manifest)


def test_control_can_run_without_releasing_mertk(tmp_path):
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({"experiments": {"pdl1_control": {"epitope": None}}}))
    check_manifest(manifest, generation=True)


def test_suite_and_child_targets_must_agree(tmp_path):
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({"experiments": {"Ig2_L270": {"epitope": "Ig2_L270"}}}))
    child = tmp_path / "campaigns/Ig2_L270/manifest.json"
    child.parent.mkdir(parents=True)
    child.write_text(json.dumps({"epitope": "Ig1_E152", "hotspots": [122, 148]}))
    with pytest.raises(ValueError, match="identities disagree"):
        check_manifest(manifest)


def test_rfd4_repair_requires_hash_matched_glycan_conditioning(tmp_path):
    import hashlib

    manifest = tmp_path / "manifest.json"
    spec = {"epitope": "Ig1_E152", "annotated_input": "input.cif"}
    manifest.write_text(
        json.dumps(
            {
                "workflow": "rfd4_native_repair_v1",
                "generation_glycan_mode": "fixed_experimental_glycan_coordinates",
                "experiments": {"parent": spec},
            }
        )
    )
    child = tmp_path / "campaigns/parent/manifest.json"
    child.parent.mkdir(parents=True)
    child.write_text(json.dumps({"epitope": "Ig1_E152", "hotspots": [122, 148]}))
    cif = tmp_path / "input.cif"
    cif.write_text("frozen input fixture")
    report = {
        "input_sha256": hashlib.sha256(cif.read_bytes()).hexdigest(),
        "pipeline_passed": True,
        "roundtrip_conditions_exact": True,
        "glycan_atoms": 520,
        "fixed_atoms": 3500,
    }
    audit = tmp_path / "input_validation.json"
    audit.write_text(json.dumps({"parent": report}))
    check_manifest(manifest, generation=True)
    cif.write_text("changed conditions")
    with pytest.raises(ValueError, match="input validation"):
        check_manifest(manifest, generation=True)
    report["input_sha256"] = hashlib.sha256(cif.read_bytes()).hexdigest()
    report["glycan_atoms"] = 0
    audit.write_text(json.dumps({"parent": report}))
    with pytest.raises(ValueError, match="input validation"):
        check_manifest(manifest, generation=True)
