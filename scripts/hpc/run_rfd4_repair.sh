#!/usr/bin/env bash
set -euo pipefail
stage=${1:?author, generate or quality}
project=${2:?project}
suite=${3:?suite}
module purge
export PYTHONNOUSERSITE=1
export PYTHONPATH="$suite/workflow"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4} MKL_NUM_THREADS=1
export HF_HOME="$project/cache/huggingface" XDG_CACHE_HOME="$project/cache" TORCH_HOME="$project/cache/torch"
export DATA_PATH="$project/cache/rfd4-data" CLUSTER=local
# PyTorch multiprocessing uses Unix-domain sockets under TMPDIR.
export TMPDIR
TMPDIR=$(mktemp -d "/tmp/rfd4-${SLURM_JOB_ID:-manual}-XXXXXX")
trap 'rm -rf -- "$TMPDIR"' EXIT
cd "$project/external/RFD4-Proteina"
source .venv/bin/activate
python - "$suite" <<'PY'
import hashlib,json,sys,subprocess
from pathlib import Path
s=Path(sys.argv[1]);r=json.loads((s/'manifest.json').read_text())
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==r['tools']['rfd4_proteina']['revision']
for name,expected in r['input_sha256'].items():
 assert hashlib.sha256((s/name).read_bytes()).hexdigest()==expected,name
print('FROZEN_RFD4_REPAIR_INPUTS_VERIFIED',flush=True)
PY
case "$stage" in
 author)
  python "$suite/workflow/rfd4_repair.py" author "$suite"
  for input in "$suite/inputs/Ig1_mpnn12.cif" "$suite/inputs/Ig1_mpnn8.cif"; do
   bash skills/rfd4-proteina-annotated-cif/scripts/inspect_cif.sh "$input"
   bash skills/rfd4-proteina-annotated-cif/scripts/check_annotated_cif.sh "$input"
  done
  ;;
 generate)
  "$project/envs/project/bin/python" "$project/src/mertk/target_policy.py" "$suite/manifest.json" --generation
  bash skills/rfd4-proteina-setup/scripts/preflight.sh
  nvidia-smi > "$suite/gpu.txt"
  for key in Ig1_mpnn12 Ig1_mpnn8; do
   rfd4-generate "inputs=$suite/inputs/$key.cif" recycle_strategy=medium/front_loaded \
    generation.args.partial_trajectory_fraction_steps=0.3 num_replicates=2 batch_size=1 \
    num_workers=2 seed=111720263 "out_dir=$suite/$key" '~metrics' compile_lengths=null \
    "hydra.run.dir=$suite/$key/hydra" prevalidate_inputs=true save_chain_of_interest=true skip_existing=false
  done
  ;;
 quality)
  "$project/envs/project/bin/python" "$suite/workflow/rfd4_repair.py" quality "$suite"
  ;;
 *) exit 2 ;;
esac
