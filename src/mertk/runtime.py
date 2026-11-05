"""Device selection, precision, synchronization, and accelerator RNG state."""

from contextlib import nullcontext

import torch


def resolve_device(requested="auto"):
    if requested == "auto":
        requested = (
            "cuda"
            if torch.cuda.is_available()
            else "mps"
            if torch.backends.mps.is_available()
            else "cpu"
        )
    device = torch.device(requested)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA is unavailable; select --device mps, cpu, or auto")
    if device.type == "mps" and not torch.backends.mps.is_available():
        raise RuntimeError(
            "MPS is unavailable to this process. Run from an unrestricted terminal on a supported Mac, or select --device cpu."
        )
    if device.type not in {"cuda", "mps", "cpu"}:
        raise ValueError("Supported devices are auto, cuda, mps, and cpu")
    return device


def resolve_precision(device, requested="auto"):
    if requested == "auto":
        requested = "bf16" if device.type == "cuda" and torch.cuda.is_bf16_supported() else "fp32"
    if requested not in {"fp32", "bf16"}:
        raise ValueError("Supported precision modes are auto, fp32, and bf16")
    if requested == "bf16":
        # ESM3 checkpoint/resume regression on MPS + bf16 produced corrupt boolean
        # masks and non-finite validation on torch 2.11. Keep this path disabled
        # until that integration passes, even though MPS supports bf16 in general.
        if device.type == "mps":
            raise RuntimeError(
                "MPS bf16 is disabled after a failed checkpoint-resume validation; use --precision fp32"
            )
        if device.type == "cpu":
            raise ValueError("Use fp32 on CPU; bf16 is enabled here only for CUDA")
        if device.type == "cuda" and not torch.cuda.is_bf16_supported():
            raise RuntimeError("This CUDA device does not support bf16; select --precision fp32")
    return requested


def autocast(device, precision="auto"):
    precision = resolve_precision(device, precision)
    return (
        torch.autocast(device.type, dtype=torch.bfloat16) if precision == "bf16" else nullcontext()
    )


def synchronize(device):
    if device.type == "mps":
        torch.mps.synchronize()
    elif device.type == "cuda":
        torch.cuda.synchronize(device)


def accelerator_rng_state(device):
    if device.type == "mps":
        return torch.mps.get_rng_state()
    if device.type == "cuda":
        return torch.cuda.get_rng_state_all()
    return None


def restore_accelerator_rng(device, state):
    if device.type == "mps":
        if state is None:
            raise ValueError("Checkpoint is missing MPS RNG state; exact resume is unavailable")
        torch.mps.set_rng_state(state)
    elif device.type == "cuda":
        if state is None:
            raise ValueError("Checkpoint is missing CUDA RNG state")
        torch.cuda.set_rng_state_all(state)
