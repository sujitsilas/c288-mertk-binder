#!/usr/bin/env bash
set -euo pipefail
stage=${1:?stage}
project=${2:?project}
suite=${3:?suite}
experiment=${4:-}
module purge
export PYTHONNOUSERSITE=1 PYTHONPATH="$suite/workflow" MPLBACKEND=Agg
export XDG_CACHE_HOME="$project/cache" HF_HOME="$project/cache/huggingface"
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export OPENBLAS_NUM_THREADS=$OMP_NUM_THREADS MKL_NUM_THREADS=$OMP_NUM_THREADS
export TMPDIR="$project/tmp/protector-${SLURM_JOB_ID:-manual}"
mkdir -p "$TMPDIR"
py="$project/envs/project/bin/python"
"$py" "$suite/workflow/bindcraft_round.py" validate "$suite" --remote
case "$stage" in
  validate|preflight)
    "$py" "$project/src/mertk/target_policy.py" "$suite/manifest.json" --generation
    "$py" "$suite/workflow/protector_round.py" "$stage" "$suite"
    ;;
  generate)
    "$py" "$project/src/mertk/target_policy.py" "$suite/manifest.json" --experiment "$experiment" --generation
    "$py" "$suite/workflow/protector_round.py" check-preflight "$suite" --experiment "$experiment"
    campaign="$suite/campaigns/$experiment"
    read -r attempts seed < <("$py" - "$campaign/manifest.json" <<'PY'
import json,sys
d=json.load(open(sys.argv[1]))['specification'];print(d['max_attempts'],d['seed'])
PY
)
    mkdir -p "$campaign/bindcraft"
    nvidia-smi > "$campaign/bindcraft/gpu.txt"
    JAX_PLATFORMS=cuda "$project/envs/bindcraft-wheels/bin/python" -u "$suite/workflow/run_bounded_bindcraft.py" \
      --root "$project/external/BindCraft" --maximum-attempts "$attempts" --seed "$seed" \
      --settings "$campaign/inputs/bindcraft/target.json" --advanced "$campaign/inputs/bindcraft/advanced.json" \
      --filters "$campaign/inputs/bindcraft/filters.json" --glycan-constraints "$campaign/inputs/glycan_constraints.json" \
      2>&1 | tee "$campaign/bindcraft/run.log"
    ;;
  quality|select|af3|report)
    "$py" "$suite/workflow/protector_round.py" "$stage" "$suite" --experiment "$experiment"
    ;;
  *) exit 2 ;;
esac
