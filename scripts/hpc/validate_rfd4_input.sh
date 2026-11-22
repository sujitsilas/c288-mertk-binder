#!/usr/bin/env bash
set -euo pipefail
project=${1:?remote project directory}
campaign=${2:?remote campaign directory}
module purge
export PYTHONNOUSERSITE=1
export CLUSTER=local
export HF_HOME="$project/cache/huggingface"
export DATA_PATH="$project/cache/rfd4-data"
export XDG_CACHE_HOME="$project/cache"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export TMPDIR="$project/tmp/input-check-${SLURM_JOB_ID:-manual}"
mkdir -p "$TMPDIR"
source "$project/external/RFD4-Proteina/.venv/bin/activate"
cd "$project/external/RFD4-Proteina"
python scripts/write_condition_spec_cifs.py "$campaign/inputs/rfd4.yaml" \
  --out-dir "$campaign/inputs/cpu_validation"
input="$campaign/inputs/cpu_validation/mertk_surface_pilot.cif"
bash skills/rfd4-proteina-annotated-cif/scripts/inspect_cif.sh "$input"
bash skills/rfd4-proteina-annotated-cif/scripts/check_annotated_cif.sh "$input"
