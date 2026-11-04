#!/usr/bin/env python3
"""Fetch a bounded set of expanded biological assemblies, plus TAM holdout seeds."""

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

from mertk.pdb import assembly_query, discover, download_assembly, download_tam_seeds, session


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=20)
    parser.add_argument("--resolution", type=float, default=3.0)
    parser.add_argument(
        "--assemblies", nargs="+", help="Explicit IDs, e.g. 1BRS-1; bypasses search filters"
    )
    parser.add_argument("--output", type=Path, default=Path("data/raw"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    client = session()
    query = assembly_query(args.resolution)
    identifiers = args.assemblies or discover(client, query, args.limit)
    records, failures = [], []
    for identifier in identifiers:
        try:
            records.append(download_assembly(client, identifier, args.output))
            print(f"Downloaded {identifier}", flush=True)
        except (OSError, ValueError) as error:
            failures.append({"assembly": identifier, "error": str(error)})
            print(f"Skipped {identifier}: {error}", flush=True)
    seeds = download_tam_seeds(client, args.output / "holdout")
    manifest = {
        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "query": None if args.assemblies else query,
        "explicit_ids": args.assemblies,
        "assemblies": records,
        "failures": failures,
        "tam_seeds": str(seeds),
    }
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    if not records:
        raise SystemExit("No assemblies downloaded; see manifest.json")


if __name__ == "__main__":
    main()
