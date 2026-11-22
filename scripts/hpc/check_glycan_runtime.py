#!/usr/bin/env python3
"""Numerical/JAX checks and optional one-step AF2 GPU integration preflight."""

import argparse
import importlib.metadata
import json
import os
from datetime import datetime, timezone
from pathlib import Path

# Bound CPU worker creation to the allocation (also permits a tiny login-node test).
if hasattr(os, "sched_getaffinity"):
    allowed_cpus = sorted(os.sched_getaffinity(0))
    os.sched_setaffinity(0, allowed_cpus[: int(os.environ.get("SLURM_CPUS_PER_TASK", "2"))])

import jax
import jax.numpy as jnp
import numpy as np

from mertk.glycan_guard import attach, losses
from mertk.provenance import file_hash


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--constraints", type=Path, required=True)
    parser.add_argument("--target", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--params", type=Path)
    parser.add_argument("--gpu-forward", action="store_true")
    parser.add_argument("--hotspot", default="A270")
    args = parser.parse_args()
    spec = json.loads(args.constraints.read_text())
    if file_hash(args.target) != spec["target_sha256"]:
        raise ValueError("Target checksum mismatch")
    xyz = np.zeros((2, 37, 3), dtype=np.float32)
    xyz[0, :3] = [[0, 1, 0], [0, 0, 0], [1, 0, 0]]
    xyz[1, :] = [0.4, 0, 3]
    mask = np.zeros((2, 37), dtype=np.float32)
    mask[:, :3] = 1
    test = {
        "target_length": 1,
        "anchor_indices": [0],
        "sequon_indices": [0],
        "glycan_local_xyz": [[0, 0, 3]],
        "glycan_atom_sites": [0],
        "heavy_atom_clearance_A": 3.5,
        "sequon_ca_exclusion_A": 15.0,
    }

    def objective(p):
        return sum(losses(p, jnp.asarray(mask), test, jnp).values())

    value, grad = jax.jit(jax.value_and_grad(objective))(jnp.asarray(xyz))
    if not np.isfinite(grad).all() or float(value) <= 0:
        raise ValueError("Invalid collision gradient")
    moved = xyz.copy()
    moved[1] -= 0.1 * np.asarray(grad)[1]
    if not float(objective(moved)) < float(value):
        raise ValueError("Gradient step does not repel the colliding binder")
    far = xyz.copy()
    far[1] += 100
    if float(objective(far)) != 0:
        raise ValueError("Distant control still penalized")
    membrane_checked = False
    if "membrane" in spec:
        p = np.zeros((4, 37, 3), dtype=np.float32)
        for i in range(3):
            p[i, :3] = np.array([[0, 1, 0], [0, 0, 0], [1, 0, 0]]) + [0, 0, i * 5]
        p[3, :3] = [30, 0, 1]
        m = np.zeros((4, 37), dtype=np.float32)
        m[:, :3] = 1
        s = dict(
            test,
            target_length=3,
            membrane={
                "outer_indices": [0],
                "inner_indices": [2],
                "origin_index": 0,
                "extracellular_margin_A": 5.0,
            },
        )
        def fn(z):
            return losses(z, jnp.asarray(m), s, jnp)["membrane_volume"]
        value_m, grad_m = jax.jit(jax.value_and_grad(fn))(jnp.asarray(p))
        moved_m = p.copy()
        moved_m[3] -= 0.1 * np.asarray(grad_m)[3]
        if not np.isfinite(grad_m).all() or not float(fn(moved_m)) < float(value_m):
            raise ValueError("Membrane gradient does not repel the binder")
        p[3, :, 2] = -20
        if float(fn(p)) != 0:
            raise ValueError("Extracellular membrane control is penalized")
        membrane_checked = True
    report = {
        "passed": True,
        "checked_at": datetime.now(timezone.utc).isoformat(),
        "constraints_sha256": file_hash(args.constraints),
        "target_sha256": file_hash(args.target),
        "jax_version": importlib.metadata.version("jax"),
        "devices": [str(d) for d in jax.devices()],
        "finite_gradient": True,
        "repulsive_gradient": True,
        "distant_control_zero": True,
        "full_forward_backward": False,
        "membrane_gradient_checked": membrane_checked,
    }
    if args.gpu_forward:
        if not any(d.platform == "gpu" for d in jax.devices()):
            raise ValueError("GPU preflight requires a GPU")
        from colabdesign import mk_afdesign_model

        model = mk_afdesign_model(
            protocol="binder",
            use_multimer=True,
            model_names=["model_1_multimer_v3"],
            num_recycles=0,
            data_dir=str(args.params),
        )
        model.prep_inputs(
            pdb_filename=str(args.target),
            chain="A",
            binder_len=65,
            hotspot=args.hotspot,
            seed=11142026,
            rm_aa="C",
            rm_target_seq=False,
            rm_target_sc=False,
        )
        attach(model, spec, args.target)
        model.design_logits(1, models=[0], num_models=1, sample_models=False, save_best=False)
        measured = {k: float(model.aux["losses"][k]) for k in spec["weights"]}
        if (
            not all(np.isfinite(v) for v in measured.values())
            or not np.isfinite(model._params["seq"]).all()
        ):
            raise ValueError("Nonfinite AF2 custom loss or updated sequence")
        report.update(full_forward_backward=True, model_losses=measured)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report), flush=True)


if __name__ == "__main__":
    main()
