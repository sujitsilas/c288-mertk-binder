#!/usr/bin/env bash
# Scheduled CPU-only dependency check. No installation or model execution.
set -euo pipefail
project=${1:?scratch project root}
module purge
export PYTHONNOUSERSITE=1
export UV_CACHE_DIR="$project/cache/uv"
export PIP_CACHE_DIR="$project/cache/pip"
export XDG_CACHE_HOME="$project/cache"
export TMPDIR="$project/tmp/bindcraft-wheel-check-${SLURM_JOB_ID:-manual}"
mkdir -p "$TMPDIR" "$project/provenance"
uv="$project/envs/alphafold3-bootstrap/bin/uv"
python="$project/envs/miniforge/envs/BindCraft/bin/python"
"$python" -c 'import sys; assert sys.version_info[:2] == (3, 10)'
date -u
timeout --signal=TERM --kill-after=15s 600s "$uv" pip install \
  --python "$python" --dry-run -r "$project/configs/bindcraft_python310.txt" \
  --find-links https://west.rosettacommons.org/pyrosetta/quarterly/release.cxx11thread.serialization
timeout --signal=TERM --kill-after=15s 300s "$uv" pip install \
  --python "$python" --dry-run --no-deps \
  'git+https://github.com/sokrypton/ColabDesign.git@8fd167defd1b0f98e29f00584ff4f3987e6935dc'
printf '%s\n' 'WHEEL_DEPENDENCY_CHECK_PASSED (no packages installed)'
date -u
