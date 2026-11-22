#!/usr/bin/env python3
"""Run pinned BindCraft with a cap on loop attempts, including rejected trajectories."""

import argparse
import importlib
import json
import os
import random
import runpy
import sys
from pathlib import Path


def attempt_gate(original_check, maximum):
    """The upstream check counts relaxed PDBs; also cap unsuccessful loop attempts."""
    if maximum < 1:
        raise ValueError("Attempt limit must be positive")
    attempts = 0

    def check(*args, **kwargs):
        nonlocal attempts
        if original_check(*args, **kwargs):
            return True
        if attempts >= maximum:
            print(f"PROJECT_ATTEMPT_CAP_REACHED {attempts}/{maximum}", flush=True)
            return True
        attempts += 1
        print(f"PROJECT_ATTEMPT {attempts}/{maximum}", flush=True)
        return False

    return check


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--maximum-attempts", type=int, required=True)
    parser.add_argument("--seed", type=int, help="Seed the upstream trajectory/length sampler")
    parser.add_argument("--settings", type=Path, required=True)
    parser.add_argument("--advanced", type=Path, required=True)
    parser.add_argument("--filters", type=Path, required=True)
    parser.add_argument("--glycan-constraints", type=Path)
    args = parser.parse_args()
    settings, advanced, filters = [
        p.resolve() for p in [args.settings, args.advanced, args.filters]
    ]
    record = json.loads(settings.read_text())
    design_path = Path(record["design_path"])
    if design_path.exists() and any(design_path.iterdir()):
        raise SystemExit("Bounded pilot requires a fresh design directory; preserve prior runs")
    if args.maximum_attempts < 1:
        parser.error("Attempt limit must be positive")
    root = args.root.resolve()
    constraints = args.glycan_constraints.resolve() if args.glycan_constraints else None
    os.chdir(root)
    sys.path.insert(0, str(root))
    functions = importlib.import_module("functions")
    if args.glycan_constraints:
        from mertk.glycan_guard import install

        install(functions, constraints)
    if args.seed is not None:
        import numpy as np

        random.seed(args.seed)
        np.random.seed(args.seed)
        print(f"PROJECT_SAMPLER_SEED {args.seed}", flush=True)
    functions.check_n_trajectories = attempt_gate(
        functions.check_n_trajectories, args.maximum_attempts
    )
    sys.argv = [
        str(root / "bindcraft.py"),
        "--settings",
        str(settings),
        "--advanced",
        str(advanced),
        "--filters",
        str(filters),
    ]
    runpy.run_path(str(root / "bindcraft.py"), run_name="__main__")
    print("PROJECT_RUN_COMPLETED", flush=True)


if __name__ == "__main__":
    main()
