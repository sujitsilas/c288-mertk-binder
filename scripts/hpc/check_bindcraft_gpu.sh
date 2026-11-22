#!/usr/bin/env bash
# Bounded runtime validation; does not generate or rank binders.
set -euo pipefail
project=${1:?scratch project root}
module purge
export PYTHONNOUSERSITE=1
unset PYTHONPATH
export XDG_CACHE_HOME="$project/cache"
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export JAX_PLATFORMS=cuda
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
cd "$project/external/BindCraft"
nvidia-smi
"$project/envs/bindcraft-wheels/bin/python" - "$project" <<'PYTHON'
import json,sys
from pathlib import Path
import jax
import jax.numpy as jnp
import functions
p=Path(sys.argv[1])
assert any(d.platform == "gpu" for d in jax.devices()), jax.devices()
a=jnp.eye(512,dtype=jnp.float32)
z=jax.jit(lambda x:x@x)(a).block_until_ready()
assert bool(jnp.allclose(z,a))
x=jnp.ones((1,16,16,4),dtype=jnp.float32)
k=jnp.ones((3,3,4,8),dtype=jnp.float32)
y=jax.jit(lambda x,k:jax.lax.conv_general_dilated(x,k,(1,1),"SAME",dimension_numbers=("NHWC","HWIO","NHWC")))(x,k).block_until_ready()
assert y.shape==(1,16,16,8) and bool(jnp.isfinite(y).all())
for i in range(1,6):
 for name in [f"params_model_{i}_ptm.npz",f"params_model_{i}_multimer_v3.npz"]:
  assert (p/"external/BindCraft/params"/name).stat().st_size>1000000,name
ready=p/"provenance/bindcraft-wheels-ready.json"
r=json.loads(ready.read_text());r.update(gpu_validated=True,gpu_devices=[str(d) for d in jax.devices()],gpu_checks=["BindCraft function imports","JIT matrix multiplication","JIT convolution","PTM/multimer parameter files present"])
ready.write_text(json.dumps(r,indent=2)+"\n")
print("BINDCRAFT_GPU_CHECK_PASSED",jax.devices())
PYTHON
