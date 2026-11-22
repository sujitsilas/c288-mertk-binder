#!/usr/bin/env bash
# Each invocation is a bounded stage; control QC failure blocks the RFD4 MERTK dependency.
set -euo pipefail
stage=${1:?validate, bindcraft, rfd4_control, rfd4_proteina or quality}
project=${2:?scratch project directory}
campaign=${3:?fresh campaign directory}
module purge
export PYTHONNOUSERSITE=1
unset PYTHONPATH
export HF_HOME="$project/cache/huggingface"
export XDG_CACHE_HOME="$project/cache"
export TORCH_HOME="$project/cache/torch"
export DATA_PATH="$project/cache/rfd4-data"
export CLUSTER=local
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export OPENBLAS_NUM_THREADS=$OMP_NUM_THREADS
export MKL_NUM_THREADS=$OMP_NUM_THREADS
export MPLBACKEND=Agg
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export TMPDIR="$project/tmp/course-$stage-${SLURM_JOB_ID:-manual}"
mkdir -p "$TMPDIR"
cd "$project"
if [[ "$stage" == bindcraft || "$stage" == rfd4_proteina ]]; then
  "$project/envs/project/bin/python" "$project/src/mertk/target_policy.py" \
    "$campaign/manifest.json" --generation
fi
"$project/envs/project/bin/python" - "$campaign" <<'PY'
import hashlib,json,sys
from pathlib import Path
campaign=Path(sys.argv[1])
record=json.loads((campaign/"manifest.json").read_text())
for relative, expected in record["input_sha256"].items():
    path=campaign/relative
    actual=hashlib.sha256(path.read_bytes()).hexdigest()
    if actual != expected:
        raise ValueError(f"Campaign input changed: {relative}")
print("CAMPAIGN_INPUT_HASHES_VERIFIED", flush=True)
PY
date -u
case "$stage" in
  validate)
    source "$project/external/RFD4-Proteina/.venv/bin/activate"
    cd "$project/external/RFD4-Proteina"
    for spec in rfd4 rfd4_control; do
      python scripts/write_condition_spec_cifs.py "$campaign/inputs/$spec.yaml" --out-dir "$campaign/inputs/$spec-annotated"
      for input in "$campaign/inputs/$spec-annotated/"*.cif; do
        bash skills/rfd4-proteina-annotated-cif/scripts/inspect_cif.sh "$input"
        bash skills/rfd4-proteina-annotated-cif/scripts/check_annotated_cif.sh "$input"
      done
    done
    python "$project/scripts/hpc/audit_course_conditions.py" \
      "$campaign/inputs/rfd4-annotated/mertk_design.cif" "$campaign/inputs/conditions_audit.json" \
      --reference "$campaign/inputs/target_reference.json"
    ;;
  bindcraft)
    "$project/envs/project/bin/python" - "$project" <<'PY'
import json,sys
from pathlib import Path
r=json.loads((Path(sys.argv[1])/"provenance/bindcraft-wheels-ready.json").read_text())
assert r["gpu_validated"] is True, "GPU smoke test not passed"
PY
    source "$project/envs/bindcraft-wheels/bin/activate"
    export JAX_PLATFORMS=cuda
    cd "$project/external/BindCraft"
    mkdir -p "$campaign/bindcraft"
    nvidia-smi > "$campaign/bindcraft/gpu.txt"
    python -c 'import jax; print(jax.devices()); assert any(d.platform == "gpu" for d in jax.devices())'
    maximum_attempts=$("$project/envs/project/bin/python" -c 'import json,sys; print(json.load(open(sys.argv[1])).get("bindcraft_max_attempts", 2))' "$campaign/manifest.json")
    python -u "$project/scripts/hpc/run_bounded_bindcraft.py" --root "$project/external/BindCraft" \
      --maximum-attempts "$maximum_attempts" --settings "$campaign/inputs/bindcraft/target.json" \
      --advanced "$campaign/inputs/bindcraft/advanced.json" --filters "$campaign/inputs/bindcraft/filters.json"
    ;;
  rfd4_control|rfd4_proteina)
    source "$project/external/RFD4-Proteina/.venv/bin/activate"
    cd "$project/external/RFD4-Proteina"
    bash skills/rfd4-proteina-setup/scripts/preflight.sh
    if [[ "$stage" == rfd4_control ]]; then
      input="$campaign/inputs/rfd4_control-annotated/pdl1_binder.cif"
      count=1; seed=101
    else
      "$project/envs/project/bin/python" - "$campaign" <<'PY'
import json,sys
from pathlib import Path
assert json.loads((Path(sys.argv[1])/"rfd4_control/quality.json").read_text())["passed"] is True
PY
      input="$campaign/inputs/rfd4-annotated/mertk_design.cif"
      count=2; seed=202
    fi
    mkdir -p "$campaign/$stage"
    nvidia-smi > "$campaign/$stage/gpu.txt"
    rfd4-generate "inputs=$input" recycle_strategy=medium/front_loaded \
      "num_replicates=$count" batch_size=1 num_workers=2 "seed=$seed" \
      "out_dir=$campaign/$stage" '~metrics' compile_lengths=null \
      "hydra.run.dir=$campaign/$stage/hydra" prevalidate_inputs=true \
      save_chain_of_interest=true skip_existing=false
    if [[ "$stage" == rfd4_control ]]; then
      cd "$project"
      "$project/envs/project/bin/python" scripts/hpc/check_course_outputs.py "$campaign" --control
    fi
    ;;
  quality)
    "$project/envs/project/bin/python" scripts/hpc/check_course_outputs.py "$campaign"
    ;;
  *) echo "Unknown stage: $stage" >&2; exit 2 ;;
esac
date -u
