import pytest
import torch

from mertk.runtime import (
    accelerator_rng_state,
    resolve_device,
    resolve_precision,
    restore_accelerator_rng,
)


def test_auto_device_prefers_cuda_then_mps_then_cpu(monkeypatch):
    monkeypatch.setattr(torch.cuda, "is_available", lambda: True)
    monkeypatch.setattr(torch.backends.mps, "is_available", lambda: True)
    assert resolve_device().type == "cuda"
    monkeypatch.setattr(torch.cuda, "is_available", lambda: False)
    assert resolve_device().type == "mps"
    monkeypatch.setattr(torch.backends.mps, "is_available", lambda: False)
    assert resolve_device().type == "cpu"
    with pytest.raises(RuntimeError, match="MPS is unavailable"):
        resolve_device("mps")


def test_precision_validation_and_conservative_defaults(monkeypatch):
    assert resolve_precision(torch.device("mps")) == "fp32"
    assert resolve_precision(torch.device("cpu")) == "fp32"
    with pytest.raises(ValueError, match="Use fp32"):
        resolve_precision(torch.device("cpu"), "bf16")
    with pytest.raises(RuntimeError, match="MPS bf16 is disabled"):
        resolve_precision(torch.device("mps"), "bf16")


def test_mps_rng_is_saved_and_restored(monkeypatch):
    expected = torch.tensor([42], dtype=torch.uint8)
    restored = []
    monkeypatch.setattr(torch.mps, "get_rng_state", lambda: expected)
    monkeypatch.setattr(torch.mps, "set_rng_state", restored.append)
    state = accelerator_rng_state(torch.device("mps"))
    restore_accelerator_rng(torch.device("mps"), state)
    assert torch.equal(restored[0], expected)
    with pytest.raises(ValueError, match="missing MPS RNG"):
        restore_accelerator_rng(torch.device("mps"), None)
