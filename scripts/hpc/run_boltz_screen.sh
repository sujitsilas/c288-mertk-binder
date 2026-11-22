#!/usr/bin/env bash
set -euo pipefail
stage=${1:?prepare or predict}
project=${2:?project directory}
suite=${3:?run directory}
module purge
module load python/3.10.20
export PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1 OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export TMPDIR
TMPDIR=$(mktemp -d "/tmp/mbz-${SLURM_JOB_ID:-local}-XXXXXX")
trap 'rm -rf -- "$TMPDIR"' EXIT
runtime="$project/envs/boltz-2.2.1-cu128"
cache="$project/cache/boltz-2-6fdef46"
cd "$suite"
sha256sum -c inputs.sha256
if [[ "$stage" == prepare ]]; then
    uv="$project/envs/rfd4_proteina-bootstrap/bin/uv"
    export UV_CACHE_DIR="$project/cache/uv" UV_LINK_MODE=hardlink
    mkdir -p "$cache"
    if [[ ! -x "$runtime/bin/python" ]]; then "$uv" venv --python "$(command -v python3)" "$runtime"; fi
    "$uv" pip install --python "$runtime/bin/python" --torch-backend cu128 'torch==2.8.0' 'boltz==2.2.1' 'setuptools<81'
    "$uv" pip check --python "$runtime/bin/python"
    "$uv" pip freeze --python "$runtime/bin/python" > runtime.txt
    "$runtime/bin/boltz" predict --help >/dev/null
    while read -r expected name; do
        if [[ ! -f "$cache/$name" ]]; then
            curl --fail --location --retry 3 --silent --show-error --output "$cache/$name.part" "https://huggingface.co/boltz-community/boltz-2/resolve/6fdef46d763fee7fbb83ca5501ccceff43b85607/$name"
            echo "$expected  $cache/$name.part" | sha256sum -c -
            mv "$cache/$name.part" "$cache/$name"
        else
            echo "$expected  $cache/$name" | sha256sum -c -
        fi
    done < weights.sha256
    if [[ ! -d "$cache/mols" ]]; then
        staging=$(mktemp -d "$cache/extract-XXXXXX")
        tar -xf "$cache/mols.tar" -C "$staging"
        mv "$staging/mols" "$cache/mols"
        rmdir "$staging"
    fi
    "$runtime/bin/python" - "$cache" <<'PY'
import json, sys
from pathlib import Path
import torch
from boltz.data.mol import load_canonicals
from boltz.data.parse.yaml import parse_yaml
mol_dir = Path(sys.argv[1]) / "mols"
canonical = load_canonicals(mol_dir)
for path in sorted(Path("inputs").glob("*.yaml")):
    record = parse_yaml(path, canonical, mol_dir, True).record
    assert len(record.chains) == 2 and all(c.msa_id == -1 for c in record.chains)
    print("Parsed", record.id, [c.num_residues for c in record.chains], flush=True)
print("CUDA build", torch.__version__, torch.version.cuda, flush=True)
assert torch.version.cuda == "12.8"
PY
elif [[ "$stage" == predict ]]; then
    nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv
    "$runtime/bin/python" -c 'import torch; assert torch.cuda.is_available() and torch.cuda.is_bf16_supported()'
    shopt -s nullglob
    inputs=(inputs/*.yaml)
    if (( ${#inputs[@]} == 0 )); then
        echo "No Boltz YAML inputs found." >&2
        exit 2
    fi
    failed=0
    prediction_complete() {
        local result=$1 name=$2
        [[ -s "$result/confidence_${name}_model_0.json" &&
            -s "$result/${name}_model_0.cif" &&
            -s "$result/pae_${name}_model_0.npz" ]]
    }
    for input in "${inputs[@]}"; do
        name=$(basename "$input" .yaml)
        result="results/boltz_results_$name/predictions/$name"
        if prediction_complete "$result" "$name"; then
            continue
        fi
        args=(predict "$input" --cache "$cache" --out_dir results --model boltz2
            --accelerator gpu --devices 1 --num_workers 0 --no_kernels --seed 101
            --recycling_steps 3 --sampling_steps 200 --diffusion_samples 1 --write_full_pae)
        # Boltz can skip a case when its prediction directory already exists.
        if [[ -d "$result" ]]; then
            args+=(--override)
        fi
        if ! "$runtime/bin/boltz" "${args[@]}" || ! prediction_complete "$result" "$name"; then
            echo "Incomplete Boltz prediction: $name" >&2
            failed=1
        fi
    done
    exit "$failed"
else
    echo "Unknown stage $stage" >&2
    exit 2
fi
