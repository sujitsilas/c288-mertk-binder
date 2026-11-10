#!/usr/bin/env python3
"""Fetch pinned upstream tools without changing the project's Python environment."""

import argparse
import json
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--tool", choices=["proteinmpnn", "bindcraft", "alphafold3", "all"], default="all"
    )
    args = parser.parse_args()
    specs = json.loads(Path("configs/external_tools.json").read_text())
    for name, spec in specs.items():
        if args.tool not in {"all", name}:
            continue
        directory = Path(spec["directory"])
        if directory.exists():
            revision = subprocess.check_output(
                ["git", "-C", str(directory), "rev-parse", "HEAD"], text=True
            ).strip()
            dirty = subprocess.check_output(
                ["git", "-C", str(directory), "status", "--porcelain", "--untracked-files=no"],
                text=True,
            )
            if revision != spec["revision"] or dirty:
                raise RuntimeError(
                    f"Existing checkout is modified or at a different revision: {directory}"
                )
        else:
            directory.parent.mkdir(parents=True, exist_ok=True)
            subprocess.run(
                ["git", "clone", "--no-checkout", spec["repository"], str(directory)], check=True
            )
            subprocess.run(
                ["git", "-C", str(directory), "checkout", "--detach", spec["revision"]], check=True
            )
        print(f"{name}: {directory} @ {spec['revision']}")
    if args.tool in {"all", "bindcraft"}:
        print(
            "BindCraft source is ready. Install its separate environment on a CUDA host using scripts/install_bindcraft.sh."
        )


if __name__ == "__main__":
    main()
