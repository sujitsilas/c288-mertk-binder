"""Releasing GPU jobs requires an unambiguous successful allocation record."""

import runpy

import pytest

validation_succeeded = runpy.run_path("scripts/hpc/submit_course_pilot.py")["validation_succeeded"]


@pytest.mark.parametrize("suffix", ["", "|"])
def test_sacct_parsable_formats(suffix):
    assert validation_succeeded(f"123|COMPLETED|0:0{suffix}\n", "123")


@pytest.mark.parametrize(
    "record",
    [
        "123|FAILED|1:0",
        "123|COMPLETED|1:0",
        "123|RUNNING|0:0",
        "124|COMPLETED|0:0",
        "123|COMPLETED|0:0\n123|COMPLETED|0:0",
        "",
    ],
)
def test_failure_missing_and_ambiguous_records_do_not_release_gpu_jobs(record):
    assert not validation_succeeded(record, "123")
