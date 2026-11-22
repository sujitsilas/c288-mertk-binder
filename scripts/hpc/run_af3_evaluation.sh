#!/usr/bin/env bash
# Submit only after candidate import, a matched assessment plan, and parameters.
set -euo pipefail
project=${1:?remote project directory}
manifest=${2:?prepared AF3 manifest}
model_dir=${3:?authorized AF3 model-parameter directory}
module purge
export PYTHONNOUSERSITE=1
export XDG_CACHE_HOME="$project/cache"
export XLA_FLAGS=--xla_gpu_enable_triton_gemm=false
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export TMPDIR="$project/tmp/af3-${SLURM_JOB_ID:-manual}"
mkdir -p "$TMPDIR"
cd "$project"
[[ -f "$manifest" && -d "$model_dir" ]]
"$project/envs/alphafold3/bin/python" -c 'import jax; print(jax.devices()); assert any(d.platform == "gpu" for d in jax.devices())'
set +e
"$project/envs/project/bin/python" scripts/run_af3.py \
  --manifest "$manifest" --af3-python "$project/envs/alphafold3/bin/python" \
  --af3-root "$project/external/alphafold3" --model-dir "$model_dir" --execute
inference_status=$?
set -e
bundle=$(dirname "$manifest")
"$project/envs/project/bin/python" scripts/score_af3.py --manifest "$manifest" \
  --predictions "$bundle/predictions" --output "$bundle/scores"
exit "$inference_status"
