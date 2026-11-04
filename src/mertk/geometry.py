"""Rotation/translation invariant local target-atom conditioning."""

import torch
from torch import nn
from torch.nn import functional as F


def local_coordinates(backbone, atoms):
    """N/CA/C define a right-handed frame; invalidate missing or collinear backbones."""
    finite = torch.isfinite(backbone).all(dim=(-1, -2))
    clean = torch.nan_to_num(backbone)
    n, ca, c = clean.unbind(-2)
    c_direction = c - ca
    e1 = F.normalize(c_direction, dim=-1, eps=1e-6)
    n_direction = n - ca
    orthogonal = n_direction - (n_direction * e1).sum(-1, keepdim=True) * e1
    valid = finite & (c_direction.norm(dim=-1) > 1e-5) & (orthogonal.norm(dim=-1) > 1e-5)
    e2 = F.normalize(orthogonal, dim=-1, eps=1e-6)
    e3 = torch.linalg.cross(e1, e2)
    frame = torch.stack((e1, e2, e3), dim=-1)
    local = torch.einsum("...ki,...ij->...kj", atoms - ca.unsqueeze(-2), frame)
    return local, valid


class GeometricAdapter(nn.Module):
    def __init__(self, d_model, hidden=128, radius=12.0, rbf_bins=16):
        super().__init__()
        self.radius = radius
        self.register_buffer("centers", torch.linspace(0, radius, rbf_bins))
        self.element = nn.Embedding(6, 8)
        self.residue = nn.Embedding(20, 16)
        self.atom_mlp = nn.Sequential(
            nn.Linear(3 + rbf_bins + 8 + 16 + 1, hidden),
            nn.SiLU(),
            nn.Linear(hidden, hidden),
            nn.SiLU(),
        )
        self.score = nn.Linear(hidden, 1)
        self.out = nn.Linear(hidden, d_model, bias=False)
        # Adding an untrained adapter initially preserves the base model exactly.
        nn.init.zeros_(self.out.weight)

    def forward(self, backbone, xyz, features, mask):
        local, valid_frame = local_coordinates(backbone.float(), xyz.float())
        distance = local.norm(dim=-1)
        mask = mask & valid_frame.unsqueeze(-1) & torch.isfinite(distance)
        local = torch.where(mask.unsqueeze(-1), local, 0)
        distance = torch.where(mask, distance, 0)
        width = self.radius / (len(self.centers) - 1)
        radial = torch.exp(-((distance.unsqueeze(-1) - self.centers) / width).square())
        atom_input = torch.cat(
            (
                local / self.radius,
                radial,
                self.element(features[..., 0]),
                self.residue(features[..., 1]),
                features[..., 2:].float(),
            ),
            dim=-1,
        )
        hidden = self.atom_mlp(atom_input.to(self.out.weight.dtype))
        # Finite softmax also for residues with no target neighbors.
        score = self.score(hidden).squeeze(-1).float().masked_fill(~mask, -1e9)
        attention = score.softmax(-1) * mask
        attention = attention / attention.sum(-1, keepdim=True).clamp_min(1e-8)
        pooled = (hidden * attention.unsqueeze(-1).to(hidden.dtype)).sum(-2)
        return self.out(pooled)
