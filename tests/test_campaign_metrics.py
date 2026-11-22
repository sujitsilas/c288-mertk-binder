import json
import runpy
from pathlib import Path

import pytest

from mertk.campaign_metrics import acceptance_interval, summarize

parse_accounting = runpy.run_path("scripts/hpc/submit_bindcraft_round.py")["parse_accounting"]
crop_clearance = runpy.run_path("scripts/hpc/prepare_bindcraft_round.py")["crop_clearance"]


def test_zero_hits_are_not_a_precise_zero_rate():
    assert acceptance_interval(0, 10)["zero_hit_upper_95_one_sided"] == pytest.approx(
        0.258866, abs=1e-6
    )
    assert acceptance_interval(0, 200)["zero_hit_upper_95_one_sided"] == pytest.approx(
        0.014867, abs=1e-6
    )
    assert acceptance_interval(0, 0)["rate"] is None


def test_ranked_copies_and_variants_do_not_multiply_independent_hits(tmp_path):
    (tmp_path / "manifest.json").write_text(json.dumps({"campaign_id": "test"}))
    root = tmp_path / "bindcraft/designs"
    (root / "Accepted/Ranked").mkdir(parents=True)
    for suffix in ["mpnn1_model1", "mpnn2_model2"]:
        (root / "Accepted" / f"binder_l80_s1_{suffix}.pdb").touch()
        (root / "Accepted/Ranked" / f"1_binder_l80_s1_{suffix}.pdb").touch()
    (root / "final_design_stats.csv").write_text(
        "Design\nbinder_l80_s1_mpnn1\nbinder_l80_s1_mpnn2\n"
    )
    (tmp_path / "bindcraft/run.log").write_text(
        "Starting trajectory: binder_l80_s1\nStarting trajectory: binder_l80_s2\n"
    )
    row = summarize(tmp_path)
    assert row["accepted_variants"] == 2
    assert row["accepted_trajectories"] == 1
    assert row["attempts_completed"] == 1
    assert row["attempts_censored"] == 1
    assert row["acceptance_on_completed_attempts"]["rate"] == 1


def test_missing_or_empty_log_is_not_reported_as_zero_attempts(tmp_path):
    (tmp_path / "manifest.json").write_text(json.dumps({"campaign_id": "test"}))
    with pytest.raises(FileNotFoundError, match="No generation log"):
        summarize(tmp_path)
    log = tmp_path / "generation.log"
    log.touch()
    with pytest.raises(ValueError, match="Generation log is empty"):
        summarize(tmp_path, log_path=log)


def test_native_log_takes_precedence_and_multiple_scheduler_logs_require_selection(tmp_path):
    (tmp_path / "manifest.json").write_text(json.dumps({"campaign_id": "test"}))
    scheduler = tmp_path / "logs"
    scheduler.mkdir()
    first = scheduler / "test_bindcraft_generate_1.log"
    first.write_text("Starting trajectory: binder_l80_s1\n")
    assert summarize(tmp_path)["attempts_started"] == 1
    (scheduler / "test_bindcraft_generate_2.log").write_text("Starting trajectory: binder_l80_s2\n")
    with pytest.raises(ValueError, match="Multiple generation logs"):
        summarize(tmp_path)
    assert summarize(tmp_path, log_path=first)["trajectory_names"] == ["binder_l80_s1"]
    native = tmp_path / "bindcraft/run.log"
    native.parent.mkdir()
    native.write_text("Starting trajectory: binder_l80_s3\n")
    row = summarize(tmp_path)
    assert row["trajectory_names"] == ["binder_l80_s3"]
    assert set(row["source_sha256"]) == {str(native)}


def test_accounting_counts_elapsed_terminal_and_full_live_reservations():
    rows = parse_accounting(
        "1|course_old|COMPLETED|0:0|3600|240|gres/gpu=1|gres/gpu=1\n"
        "2|round_new|RUNNING|0:0|600|120|gres/gpu=1|gres/gpu=1\n"
        "3|round_next|PENDING|0:0|0|60||gres/gpu:a100=1\n"
        "4|other_work|COMPLETED|0:0|9999|240|gres/gpu=1|gres/gpu=1\n"
        "5|course_cancelled|CANCELLED by 123|0:0|0|240||gres/gpu=1\n"
    )
    assert sum(r["spent_or_reserved_gpu_seconds"] for r in rows) == 14400
    assert sum(r["elapsed_gpu_seconds"] for r in rows) == 4200


def test_accounting_counts_each_gpu_array_task():
    rows = parse_accounting(
        "42_0|mertk_test|COMPLETED|0:0|300|10|gres/gpu=1|gres/gpu=1\n"
        "42_1|mertk_test|PENDING|0:0|0|10||gres/gpu=1\n"
        "43_0|mertk_cpu|PENDING|0:0|0|10||cpu=4\n"
    )
    assert [row["job_id"] for row in rows] == ["42_0", "42_1"]
    assert sum(r["spent_or_reserved_gpu_seconds"] for r in rows) == 900
    with pytest.raises(ValueError, match="Cannot account"):
        parse_accounting("42_[0-1]|mertk_test|PENDING|0:0|0|10||gres/gpu=1\n")


def test_native_crops_pass_but_original_patch_a_crop_severs_context():
    from Bio.PDB import PDBParser

    reference = Path("data/references/37NT/target.pdb")
    if not reference.is_file():
        pytest.skip("37NT reference is missing. Run scripts/fetch_reference.py to install it.")
    chain = PDBParser(QUIET=True).get_structure("m", reference)[0]["A"]
    assert crop_clearance(chain, [93, 285], [122, 147, 148])["nearest_excluded_heavy_atom_A"] > 70
    assert crop_clearance(chain, [286, 484], [469, 471, 472])["nearest_excluded_heavy_atom_A"] > 40
    with pytest.raises(ValueError, match="within 20 A"):
        crop_clearance(chain, [286, 484], [365, 367, 369])
