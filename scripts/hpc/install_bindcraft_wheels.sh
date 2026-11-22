#!/usr/bin/env bash
# Scheduled CPU installation using upstream-compatible wheels; all assets in scratch.
set -euo pipefail
project=${1:?scratch project root}
module purge
export PYTHONNOUSERSITE=1
unset PYTHONPATH
export UV_CACHE_DIR="$project/cache/uv"
export PIP_CACHE_DIR="$project/cache/pip"
export XDG_CACHE_HOME="$project/cache"
export TMPDIR="$project/tmp/bindcraft-wheels-${SLURM_JOB_ID:-manual}"
export OPENMM_CPU_THREADS=${SLURM_CPUS_PER_TASK:-8}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
mkdir -p "$TMPDIR" "$project/provenance"
uv=${UV_EXECUTABLE:-"$project/envs/alphafold3-bootstrap/bin/uv"}
base_python="$project/envs/miniforge/envs/BindCraft/bin/python"
env_dir="$project/envs/bindcraft-wheels"
"$base_python" -c 'import sys; assert sys.version_info[:2] == (3, 10)'
if [[ ! -x "$env_dir/bin/python" ]]; then "$base_python" -m venv "$env_dir"; fi
cd "$project/external/BindCraft"
[[ $(git rev-parse HEAD) == 7713aa0d0d351e4117a8befeb8541f3a8ebd3368 ]]
date -u
timeout --signal=TERM --kill-after=30s 1800s "$uv" pip install   --python "$env_dir/bin/python" -r "$project/configs/bindcraft_python310.txt"   --find-links https://west.rosettacommons.org/pyrosetta/quarterly/release.cxx11thread.serialization
timeout --signal=TERM --kill-after=30s 600s "$uv" pip install   --python "$env_dir/bin/python" --no-deps   'git+https://github.com/sokrypton/ColabDesign.git@8fd167defd1b0f98e29f00584ff4f3987e6935dc'
"$uv" pip check --python "$env_dir/bin/python"
JAX_PLATFORMS=cpu "$env_dir/bin/python" - <<'PYTHON'
import jax, numpy, colabdesign, pyrosetta, pdbfixer, openmm
assert jax.__version__ == "0.4.35"
assert numpy.__version__ == "1.26.4"
print("CPU import check passed", jax.__version__, numpy.__version__)
PYTHON
chmod +x functions/dssp functions/DAlphaBall.gcc
# The AF2 parameter files used by the pinned upstream installer.
mkdir -p params
if [[ ! -f params/params_model_5_ptm.npz ]]; then
  url=https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar
  curl -fL --retry 3 "$url" -o "$project/cache/alphafold_params_2022-12-06.tar"
  sha256sum "$project/cache/alphafold_params_2022-12-06.tar" > "$project/provenance/bindcraft-af2-archive.sha256"
  tar tf "$project/cache/alphafold_params_2022-12-06.tar" > "$project/provenance/bindcraft-af2-archive-files.txt"
  tar xf "$project/cache/alphafold_params_2022-12-06.tar" -C params
fi
[[ -f params/params_model_5_ptm.npz ]]
"$uv" pip freeze --python "$env_dir/bin/python" > "$project/provenance/bindcraft-wheels-pip.txt"
"$env_dir/bin/python" - "$project" <<'PYTHON'
import json,sys
from pathlib import Path
p=Path(sys.argv[1]);(p/"provenance/bindcraft-wheels-ready.json").write_text(json.dumps({"environment":str(p/"envs/bindcraft-wheels"),"bindcraft_revision":"7713aa0d0d351e4117a8befeb8541f3a8ebd3368","colabdesign_revision":"8fd167defd1b0f98e29f00584ff4f3987e6935dc","cpu_imports_passed":True,"gpu_validated":False},indent=2)+"\n")
PYTHON
date -u
