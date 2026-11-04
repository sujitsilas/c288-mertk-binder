import numpy as np
import torch

from mertk.data import Chain, prepare_pair
from mertk.dataset import synthetic_example
from mertk.geometry import GeometricAdapter, local_coordinates


def inputs():
    example = synthetic_example()
    return [
        torch.from_numpy(example[key]).unsqueeze(0)
        for key in ("binder_backbone", "neighbor_xyz", "neighbor_features", "neighbor_mask")
    ]


def test_local_frame_and_adapter_are_rigid_motion_invariant():
    torch.manual_seed(12)
    backbone, xyz, features, mask = inputs()
    rotation, _ = torch.linalg.qr(torch.randn(3, 3))
    rotation[:, 0] *= torch.linalg.det(rotation)
    translation = torch.tensor([11.0, -4.0, 9.0])
    transformed = (backbone @ rotation + translation, xyz @ rotation + translation)
    local, valid = local_coordinates(backbone, xyz)
    moved, moved_valid = local_coordinates(*transformed)
    torch.testing.assert_close(local, moved, atol=2e-5, rtol=2e-5)
    assert torch.equal(valid, moved_valid)
    adapter = GeometricAdapter(24)
    torch.nn.init.normal_(adapter.out.weight)
    torch.testing.assert_close(
        adapter(backbone, xyz, features, mask),
        adapter(*transformed, features, mask),
        atol=2e-5,
        rtol=2e-5,
    )


def test_empty_neighbors_and_invalid_backbones_produce_zero_finite_output():
    backbone, xyz, features, mask = inputs()
    backbone[:, 0] = float("nan")
    backbone[:, 1] = 0  # Degenerate frame.
    mask[:, 2] = False
    adapter = GeometricAdapter(24)
    torch.nn.init.normal_(adapter.out.weight)
    result = adapter(backbone, xyz, features, mask)
    assert torch.isfinite(result).all()
    assert (result[:, :3] == 0).all()
    result.sum().backward()
    assert all(torch.isfinite(p.grad).all() for p in adapter.parameters() if p.grad is not None)


def test_binder_sidechains_cannot_change_adapter_inputs_or_interface_mask():
    fixture = synthetic_example()
    binder = Chain(
        "A",
        str(fixture["binder_sequence"]),
        fixture["binder_backbone"],
        np.random.default_rng(0).normal(size=(50, 3)),
        np.zeros((50, 3), dtype=int),
    )
    target = Chain(
        "B",
        str(fixture["target_sequence"]),
        fixture["target_backbone"],
        fixture["neighbor_xyz"].reshape(-1, 3),
        fixture["neighbor_features"].reshape(-1, 3),
    )
    before = prepare_pair(binder, target)
    binder.atoms[:] = 10000
    binder.atom_features[:] = 19
    after = prepare_pair(binder, target)
    for key in before:
        np.testing.assert_equal(before[key], after[key])
