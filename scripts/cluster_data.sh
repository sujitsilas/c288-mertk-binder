#!/usr/bin/env bash
set -euo pipefail

# MMseqs2 is an external executable, separate from Python package installation.
command -v mmseqs >/dev/null || {
    echo 'Install MMseqs2 first: https://github.com/soedinglab/MMseqs2#installation' >&2
    exit 1
}

data_dir="${1:-data/processed}"
mkdir -p "$data_dir/clustering"
mmseqs version > "$data_dir/clustering/version.txt"
mmseqs easy-cluster "$data_dir/chains.fasta" "$data_dir/clustering/families" \
    "$data_dir/clustering/tmp" --min-seq-id 0.3 -c 0.8 --cov-mode 0 \
    --cluster-mode 1 -s 7.5
python scripts/split_data.py --examples "$data_dir/examples.jsonl" \
    --clusters "$data_dir/clustering/families_cluster.tsv" --output "$data_dir/splits.jsonl"
