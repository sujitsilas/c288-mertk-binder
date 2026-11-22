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
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export OPENBLAS_NUM_THREADS=$OMP_NUM_THREADS MKL_NUM_THREADS=$OMP_NUM_THREADS
export TMPDIR="$project/tmp/ig-followup-${SLURM_JOB_ID:-manual}"
mkdir -p "$TMPDIR"
py="$project/envs/project/bin/python"
"$py" "$suite/workflow/bindcraft_round.py" validate "$suite" --remote
case "$stage" in
  validate)
    JAX_PLATFORMS=cpu "$project/envs/bindcraft-wheels/bin/python" "$suite/workflow/check_glycan_runtime.py" \
      --target "$suite/campaigns/Ig1_E152/inputs/target.pdb" \
      --constraints "$suite/campaigns/Ig1_E152/inputs/glycan_constraints.json" --output "$suite/cpu_preflight.json"
    ;;
  glycan-preflight)
    JAX_PLATFORMS=cuda "$project/envs/bindcraft-wheels/bin/python" "$suite/workflow/check_glycan_runtime.py" \
      --target "$suite/campaigns/Ig1_E152/inputs/target.pdb" \
      --constraints "$suite/campaigns/Ig1_E152/inputs/glycan_constraints.json" \
      --params "$project/external/BindCraft" --gpu-forward --output "$suite/gpu_preflight.json"
    ;;
  generate)
    "$py" "$project/src/mertk/target_policy.py" "$suite/manifest.json" --experiment "$experiment" --generation
    "$py" - "$suite" "$experiment" <<'PY'
import json,sys,hashlib
from pathlib import Path
p=Path(sys.argv[1]);q=json.loads((p/'gpu_preflight.json').read_text())
c=p/'campaigns'/sys.argv[2]/'inputs/glycan_constraints.json'
if not q.get('passed') or not q.get('full_forward_backward') or q['constraints_sha256'] != hashlib.sha256(c.read_bytes()).hexdigest():
    raise SystemExit('Validated GPU glycan preflight required')
if json.loads((p/'workflow/positive_control_evidence.json').read_text()).get('passed') is not True:
    raise SystemExit('No positive-control evidence')
PY
    campaign="$suite/campaigns/$experiment"
    read -r maximum_attempts sampler_seed < <("$py" - "$campaign/manifest.json" <<'PY'
import json,sys
d=json.load(open(sys.argv[1]))['specification'];print(d['max_attempts'],d['seed'])
PY
)
    mkdir -p "$campaign/bindcraft"
    nvidia-smi > "$campaign/bindcraft/gpu.txt"
    JAX_PLATFORMS=cuda "$project/envs/bindcraft-wheels/bin/python" -u "$suite/workflow/run_bounded_bindcraft.py" \
      --root "$project/external/BindCraft" --maximum-attempts "$maximum_attempts" --seed "$sampler_seed" \
      --settings "$campaign/inputs/bindcraft/target.json" --advanced "$campaign/inputs/bindcraft/advanced.json" \
      --filters "$campaign/inputs/bindcraft/filters.json" --glycan-constraints "$campaign/inputs/glycan_constraints.json" \
      2>&1 | tee "$campaign/bindcraft/run.log"
    ;;
  quality) "$py" "$suite/workflow/bindcraft_round.py" report "$suite" --experiment "$experiment" ;;
  protector-triage) "$py" "$suite/workflow/ig_followup.py" prepare "$suite" ;;
  af3)
    export XLA_FLAGS=--xla_gpu_enable_triton_gemm=false
    "$project/envs/alphafold3/bin/python" -c 'import jax; print(jax.devices()); assert any(d.platform == "gpu" for d in jax.devices())'
    "$py" "$suite/workflow/ig_followup.py" af3 "$suite" --case "$experiment"
    ;;
  report) "$py" "$suite/workflow/ig_followup.py" report "$suite" ;;
  *) exit 2 ;;
esac
