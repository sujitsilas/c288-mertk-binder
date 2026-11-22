"""Failed upstream attempts must consume the pilot's finite attempt budget."""

import runpy
from pathlib import Path

import pytest

gate = runpy.run_path("scripts/hpc/run_bounded_bindcraft.py")["attempt_gate"]


def test_rejected_attempts_still_reach_limit():
    check = gate(lambda *args: False, 2)
    assert check({}, {}) is False
    assert check({}, {}) is False
    assert check({}, {}) is True
    assert check({}, {}) is True


def test_upstream_stopping_condition_is_preserved():
    check = gate(lambda *args: True, 2)
    assert check({}, {}) is True


def test_attempt_count_is_independent_of_upstream_relaxed_count():
    relaxed = []
    check = gate(lambda *args: len(relaxed) >= 2, 2)
    assert check({}, {}) is False
    relaxed.append("one_completed_trajectory")
    assert check({}, {}) is False
    assert check({}, {}) is True


def test_positive_limit_required():
    with pytest.raises(ValueError, match="positive"):
        gate(lambda: False, 0)


def test_launcher_applies_gate_to_upstream_import(tmp_path, monkeypatch):
    import json
    import sys
    import types

    functions = types.ModuleType("functions")
    functions.check_n_trajectories = lambda *args: False
    monkeypatch.setitem(sys.modules, "functions", functions)
    root = tmp_path / "engine"
    root.mkdir()
    (root / "bindcraft.py").write_text(
        "from functions import check_n_trajectories\n"
        "from pathlib import Path\n"
        "n = 0\n"
        "while not check_n_trajectories({}, {}):\n"
        "    n += 1\n"
        "Path('attempts.txt').write_text(str(n))\n"
    )
    settings, advanced, filters = [
        tmp_path / name for name in ["target.json", "advanced.json", "filters.json"]
    ]
    settings.write_text(json.dumps({"design_path": str(tmp_path / "fresh_designs")}))
    for path in [advanced, filters]:
        path.write_text("{}")
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_bounded_bindcraft.py",
            "--root",
            str(root),
            "--maximum-attempts",
            "2",
            "--settings",
            str(settings),
            "--advanced",
            str(advanced),
            "--filters",
            str(filters),
        ],
    )
    monkeypatch.syspath_prepend(str(root))
    runpy.run_path(
        str(Path(__file__).resolve().parents[1] / "scripts/hpc/run_bounded_bindcraft.py"),
        run_name="__main__",
    )
    assert (root / "attempts.txt").read_text() == "2"
