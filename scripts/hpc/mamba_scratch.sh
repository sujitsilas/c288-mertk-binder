#!/usr/bin/env bash
# Keep the upstream installer independent of the user's unrelated Conda channels.
set -euo pipefail
base=${MAMBA_ROOT_PREFIX:?private Miniforge root is required}
subcommand=${1:?mamba subcommand}
shift
case "$subcommand" in
  create)
    if [[ -x "$base/envs/BindCraft/bin/python" ]]; then
      "$base/envs/BindCraft/bin/python" -c 'import sys; assert sys.version_info[:2] == (3, 10)'
      echo "Reusing the existing Python 3.10 environment at $base/envs/BindCraft"
      exit 0
    fi
    exec "$base/bin/mamba" create --override-channels --strict-channel-priority \
      -c conda-forge -c nvidia "$@"
    ;;
  install)
    # Narrow the upstream ranges: the unconstrained solve timed out after four
    # hours. conda-forge's Python 3.10 JAX 0.4.35 GPU build requires CUDA >=12.6.
    # All these versions satisfy the pinned BindCraft installer's requirements.
    # Fail with a useful log before consuming the entire setup allocation again.
    constraints=("jax=0.4.35" "jaxlib=0.4.35=cuda126py310h5e1a0f3_200"
      "cuda-version=12.6" "cuda-nvcc=12.6" "cudnn=9.3"
      "numpy=1.26.4" "dm-haiku=0.0.13" "flax=0.9.0" "optax=0.2.4" "chex=0.1.87")
    echo "Resolving the bounded BindCraft JAX 0.4.35 / CUDA 12.6 environment"
    exec timeout --signal=TERM --kill-after=30s 1800s \
      "$base/bin/mamba" install --override-channels --strict-channel-priority \
      -c conda-forge -c nvidia "${constraints[@]}" "$@"
    ;;
  *) exec "$base/bin/mamba" "$subcommand" "$@" ;;
esac
