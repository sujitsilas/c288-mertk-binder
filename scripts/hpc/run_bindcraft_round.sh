#!/usr/bin/env bash
set -euo pipefail
stage=${1:?validate, generate, control-gate or report}
project=${2:?scratch project}
suite=${3:?frozen suite}
experiment=${4:-}
module purge
export PYTHONNOUSERSITE=1
export PYTHONPATH="$suite/workflow"
export XDG_CACHE_HOME="$project/cache"
export HF_HOME="$project/cache/huggingface"
export MPLBACKEND=Agg
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export OPENBLAS_NUM_THREADS=$OMP_NUM_THREADS
export MKL_NUM_THREADS=$OMP_NUM_THREADS
export TMPDIR="$project/tmp/bindcraft-round-${SLURM_JOB_ID:-manual}"
mkdir -p "$TMPDIR"
runner="$suite/workflow/bindcraft_round.py"
python_project="$project/envs/project/bin/python"
if [[ "$stage" == generate ]]; then
  "$python_project" "$project/src/mertk/target_policy.py" "$suite/manifest.json" \
    --experiment "$experiment" --generation
fi
"$python_project" "$runner" validate "$suite" --remote
case "$stage" in
  validate) ;;
  control-gate) "$python_project" "$runner" control-gate "$suite" ;;
  report) "$python_project" "$runner" report "$suite" --experiment "$experiment" ;;
  generate)
    campaign="$suite/campaigns/$experiment"
    if [[ "$experiment" != pdl1_control ]]; then
      "$python_project" - "$suite/campaigns/pdl1_control/control_gate.json" <<'PY'
import json,sys
if json.load(open(sys.argv[1])).get("passed") is not True:
    raise SystemExit("Positive-control gate did not pass")
PY
    fi
    read -r maximum_attempts sampler_seed < <("$python_project" - "$campaign/manifest.json" <<'PY'
import json,sys
d=json.load(open(sys.argv[1]))["specification"]
print(d["max_attempts"],d["seed"])
PY
)
    source "$project/envs/bindcraft-wheels/bin/activate"
    export JAX_PLATFORMS=cuda
    mkdir -p "$campaign/bindcraft"
    nvidia-smi > "$campaign/bindcraft/gpu.txt"
    python -u "$suite/workflow/run_bounded_bindcraft.py" --root "$project/external/BindCraft" \
      --maximum-attempts "$maximum_attempts" --seed "$sampler_seed" \
      --settings "$campaign/inputs/bindcraft/target.json" \
      --advanced "$campaign/inputs/bindcraft/advanced.json" \
      --filters "$campaign/inputs/bindcraft/filters.json" 2>&1 | tee "$campaign/bindcraft/run.log"
    ;;
  *) echo "Unknown stage: $stage" >&2; exit 2 ;;
esac
