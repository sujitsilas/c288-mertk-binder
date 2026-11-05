#!/usr/bin/env bash
set -euo pipefail
if [[ "$(uname -s)" != Linux ]]; then
    echo 'The pinned BindCraft environment requires a Linux CUDA host.' >&2
    exit 1
fi
command -v nvidia-smi >/dev/null || { echo 'An NVIDIA CUDA driver is required.' >&2; exit 1; }
command -v conda >/dev/null || { echo 'Conda must be installed first.' >&2; exit 1; }
repo_root="$(cd "$(dirname "$0")/.." && pwd)"
bindcraft_dir="$repo_root/external/BindCraft"
expected_revision=7713aa0d0d351e4117a8befeb8541f3a8ebd3368
[[ "$(git -C "$bindcraft_dir" rev-parse HEAD)" == "$expected_revision" ]] || {
    echo 'Run scripts/setup_external.py --tool bindcraft first.' >&2; exit 1;
}
cd "$bindcraft_dir"
# Upstream installer creates its own conda environment and downloads AF2 weights.
# CUDA version can be supplied explicitly, e.g. --cuda 12.4 --pkg_manager conda.
bash install_bindcraft.sh "$@"
conda list -n BindCraft --explicit > environment-resolved.txt
conda run -n BindCraft python -m pip freeze > pip-resolved.txt
