import importlib.util
import json
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts/hpc"


def module(name):
    sys.path.insert(0, str(SCRIPTS))
    try:
        spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
        loaded = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(loaded)
        return loaded
    finally:
        sys.path.pop(0)


def test_budget_and_dependency_gates():
    submit = module("submit_ig_followup")
    record = json.loads(Path("configs/bindcraft_round.json").read_text())
    rows = submit.job_specs(record)
    assert sum(submit.seconds(r[-1]) for r in rows if r[-2]) / 3600 == pytest.approx(6 + 1 / 12)
    named = {r[0]: r for r in rows}
    for name in ("Ig1_E152_generate", "Ig2_L270_generate"):
        assert named[name][3:5] == (["glycan_preflight"], "afterok")
    assert named["af3_parents"][3:5] == (["protector_triage"], "afterok")
    assert named["af3_tulp2"][3:5] == (["input_check"], "afterok")


def test_no_parents_cannot_become_protectors_or_af3_inputs(tmp_path):
    follow = module("ig_followup")
    (tmp_path / "manifest.json").write_text(
        json.dumps({"experiments": {"Ig2_L270": {}}, "remote_project": "/unused"})
    )
    (tmp_path / "comparisons").mkdir()
    (tmp_path / "comparisons/manifest.json").write_text(
        json.dumps({"cases": [], "protocol": {}, "target": {}})
    )
    follow.prepare_candidates(tmp_path)
    record = json.loads((tmp_path / "protector/triage.json").read_text())
    assert record["parents"] == []
    assert record["cleavage_protection_predicted"] is False
    follow.run_af3(tmp_path, "parents")
    assert (tmp_path / "assessment/execution_skipped.json").exists()
    assert not (tmp_path / "assessment/predictions").exists()


def test_missing_control_outputs_remain_missing(tmp_path):
    follow = module("ig_followup")
    bundle = tmp_path / "comparisons"
    bundle.mkdir()
    record = {
        "protocol": {},
        "target": {},
        "cases": [{"id": "mertk_mcherry", "input": "absent.json", "input_sha256": "unavailable"}],
    }
    (bundle / "manifest.json").write_text(json.dumps(record))
    follow.report_af3(tmp_path)
    result = json.loads((bundle / "results_manifest.json").read_text())
    assert len(result["errors"]) == 1
    assert result["cases"][0]["predictions"] == []
    assert json.loads((bundle / "manifest.json").read_text()) == record
