#!/usr/bin/env bash
set -euo pipefail
phase=${1:?parents or repairs}
project=${2:?project}
suite=${3:?suite}
repair=${4:?repair}
parents=(Ig1_mpnn12 Ig1_mpnn8)
module purge
export PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1
export PYTHONPATH="$repair/workflow"
export JAX_PLATFORMS=cpu XLA_PYTHON_CLIENT_PREALLOCATE=false
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8} OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLBACKEND=Agg
export TMPDIR="$project/tmp/refold-${SLURM_JOB_ID:?}"
mkdir -p "$TMPDIR"
exec "$project/envs/bindcraft-wheels/bin/python" "$suite/workflow/score_ig1_refolds.py" \
  "$phase" "${parents[${SLURM_ARRAY_TASK_ID:?}]}" "$project" "$suite" "$repair"
