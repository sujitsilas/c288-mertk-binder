#!/usr/bin/env python3
import argparse
import json
from collections import Counter
from pathlib import Path

from mertk.splits import assign_splits

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Split MMseqs2 family components without chain overlap"
    )
    parser.add_argument("--examples", type=Path, default=Path("data/processed/examples.jsonl"))
    parser.add_argument("--clusters", type=Path, required=True, help="MMseqs2 *_cluster.tsv")
    parser.add_argument("--output", type=Path, default=Path("data/processed/splits.jsonl"))
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    clusters = {}
    for line in args.clusters.read_text().splitlines():
        representative, member = line.split("\t")
        if member in clusters and clusters[member] != representative:
            raise ValueError(f"Conflicting clusters for {member}")
        clusters[member] = representative
    rows = [json.loads(line) for line in args.examples.read_text().splitlines() if line.strip()]
    rows = assign_splits(rows, clusters, seed=args.seed)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("".join(json.dumps(row) + "\n" for row in rows))
    print(dict(Counter(row["split"] for row in rows)))
