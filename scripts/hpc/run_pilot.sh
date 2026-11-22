#!/usr/bin/env bash
# Each invocation is one bounded GPU pilot; no refolding oracle is substituted.
set -euo pipefail
tool=${1:?bindcraft or rfd4_proteina}
project=${2:?remote project directory}
campaign=${3:?campaign directory}
"$project/envs/project/bin/python" "$project/src/mertk/target_policy.py" \
  "$campaign/manifest.json" --generation
module purge
export PYTHONNOUSERSITE=1
export HF_HOME="$project/cache/huggingface"
export XDG_CACHE_HOME="$project/cache"
export TORCH_HOME="$project/cache/torch"
export DATA_PATH="$project/cache/rfd4-data"
export CLUSTER=local
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export MPLBACKEND=Agg
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export TMPDIR="$project/tmp/pilot-$tool-${SLURM_JOB_ID:-manual}"
mkdir -p "$campaign/$tool" "$TMPDIR"
nvidia-smi > "$campaign/$tool/gpu.txt"
date -u
case "$tool" in
  bindcraft)
    [[ -f "$project/provenance/bindcraft-wheels-ready.json" ]]
    source "$project/envs/bindcraft-wheels/bin/activate"
    cd "$project/external/BindCraft"
    python -c 'import jax; print(jax.devices()); assert any(d.platform == "gpu" for d in jax.devices())'
    python -u bindcraft.py --settings "$campaign/inputs/bindcraft/target.json" \
      --advanced "$campaign/inputs/bindcraft/advanced.json" \
      --filters "$campaign/inputs/bindcraft/filters.json"
    ;;
  rfd4_proteina)
    source "$project/external/RFD4-Proteina/.venv/bin/activate"
    cd "$project/external/RFD4-Proteina"
    bash skills/rfd4-proteina-setup/scripts/preflight.sh
    # The longest pilot input is 392 target + 80 binder residues, batch size 1.
    # Preserve sequence and internal backbone geometry; no reviewed hotspot yet.
    python scripts/write_condition_spec_cifs.py "$campaign/inputs/rfd4.yaml" \
      --out-dir "$campaign/inputs/annotated"
    input="$campaign/inputs/annotated/mertk_surface_pilot.cif"
    bash skills/rfd4-proteina-annotated-cif/scripts/inspect_cif.sh "$input"
    bash skills/rfd4-proteina-annotated-cif/scripts/check_annotated_cif.sh "$input"
    rfd4-generate "inputs=$input" recycle_strategy=medium/front_loaded \
      num_replicates=2 batch_size=1 num_workers=4 seed=101 \
      "out_dir=$campaign/rfd4_proteina" '~metrics' compile_lengths=null \
      prevalidate_inputs=true save_chain_of_interest=true skip_existing=false
    ;;
  *) echo "Unknown runtime: $tool" >&2; exit 2 ;;
esac
date -u
