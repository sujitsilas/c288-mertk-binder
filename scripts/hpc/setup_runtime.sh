#!/usr/bin/env bash
# Run in a scheduled CPU allocation, never on the login node.
set -euo pipefail
tool=${1:?bindcraft, rfd4_proteina or alphafold3}
project=${2:?remote project directory}
module purge
module load python/3.12.13
export UV_CACHE_DIR="$project/cache/uv"
export HF_HOME="$project/cache/huggingface"
export CONDA_PKGS_DIRS="$project/cache/conda"
export XDG_CACHE_HOME="$project/cache"
export PIP_CACHE_DIR="$project/cache/pip"
export TORCH_HOME="$project/cache/torch"
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export TMPDIR="$project/tmp/setup-$tool-${SLURM_JOB_ID:-manual}"
mkdir -p "$project/envs" "$project/cache" "$project/provenance" "$TMPDIR"
bootstrap_uv() {
    local bootstrap="$project/envs/$tool-bootstrap"
    python3.12 -m venv "$bootstrap"
    "$bootstrap/bin/python" -m pip install 'pip @ https://files.pythonhosted.org/packages/3f/06/3e4420302363fbb71be3c7f35ce042f54bfca330e3e66b3a2e3c4cca579f/pip-26.2.2-py3-none-any.whl#sha256=fd7f4dfa8190e014013b9f08f3e973f4e4c26ffd179a204895533a7254813a70'
    "$bootstrap/bin/python" -m pip install 'uv==0.13.6'
    uv="$bootstrap/bin/uv"
}
date -u
case "$tool" in
  bindcraft)
    # A private base lets the upstream installer use its expected envs/BindCraft
    # location without writing to Hoffman2's shared Miniforge installation.
    base="$project/envs/miniforge"
    if [[ ! -x "$base/bin/conda" ]]; then
      installer="$project/cache/Miniforge3-26.5.3-0-Linux-x86_64.sh"
      url=https://github.com/conda-forge/miniforge/releases/download/26.5.3-0/Miniforge3-26.5.3-0-Linux-x86_64.sh
      curl -fL --retry 3 "$url" -o "$installer"
      curl -fL --retry 3 "$url.sha256" -o "$installer.sha256"
      (cd "$project/cache" && sha256sum -c "$(basename "$installer.sha256")")
      bash "$installer" -b -p "$base"
    fi
    export PATH="$base/bin:$base/condabin:$PATH"
    # Mamba can infer the wrong root when Miniforge itself is inside envs/.
    export MAMBA_ROOT_PREFIX="$base"
    export CONDA_ENVS_PATH="$base/envs"
    cd "$project/external/BindCraft"
    [[ $(git rev-parse HEAD) == 7713aa0d0d351e4117a8befeb8541f3a8ebd3368 ]]
    # Keep the small Python-only Conda environment; the full Conda solve has
    # repeatedly stalled. Install upstream-compatible CUDA wheels in isolation.
    if [[ ! -x "$base/envs/BindCraft/bin/python" ]]; then
      timeout --signal=TERM --kill-after=30s 900s \
        "$project/scripts/hpc/mamba_scratch.sh" create --name BindCraft python=3.10 -y
    fi
    bootstrap_uv
    UV_EXECUTABLE="$uv" bash "$project/scripts/hpc/install_bindcraft_wheels.sh" "$project"
    conda list -n BindCraft --explicit > "$project/provenance/bindcraft-conda.txt"
    ;;
  rfd4_proteina)
    cd "$project/external/RFD4-Proteina"
    [[ $(git rev-parse HEAD) == 68058693cb88762a8aaf806571e13b291ce6318a ]]
    bootstrap_uv
    "$uv" venv --python python3.12 .venv
    source .venv/bin/activate
    "$uv" pip install --no-sources --torch-backend cu126 -e '.[all-cuda12]'
    "$uv" pip check
    export DATA_PATH="$project/cache/rfd4-data"
    export CLUSTER=local
    [[ -e .env.local ]] || cp .env.local.example .env.local
    bash scripts/download_assets.sh "$DATA_PATH" --skip-of3 --skip-benchmarks
    "$uv" pip freeze > "$project/provenance/rfd4-pip.txt"
    ;;
  alphafold3)
    cd "$project/external/alphafold3"
    [[ $(git rev-parse HEAD) == 28c2d290559fe07387a12f471ad49f0d36206327 ]]
    bootstrap_uv
    export UV_PROJECT_ENVIRONMENT="$project/envs/alphafold3"
    export UV_LINK_MODE=copy
    command -v g++
    "$uv" sync --frozen --no-dev --no-editable
    "$UV_PROJECT_ENVIRONMENT/bin/build_data"
    "$uv" pip freeze --python "$UV_PROJECT_ENVIRONMENT/bin/python" > "$project/provenance/af3-pip.txt"
    # The project runner has its own lightweight Python environment.
    "$uv" venv --python python3.12 "$project/envs/project"
    "$uv" pip install --python "$project/envs/project/bin/python" -e "$project"
    ;;
  *) echo "Unknown runtime: $tool" >&2; exit 2 ;;
esac
date -u
