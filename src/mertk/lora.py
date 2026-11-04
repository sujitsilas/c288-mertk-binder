"""LoRA for ESM3's fused QKV and attention output projections."""

import math

import torch
from torch import nn


class LoRALinear(nn.Module):
    def __init__(self, base, rank=8, alpha=16, dropout=0.05):
        super().__init__()
        if rank < 1:
            raise ValueError("LoRA rank must be positive")
        self.base = base.requires_grad_(False)
        self.scale = alpha / rank
        self.dropout = nn.Dropout(dropout)
        self.lora_a = nn.Parameter(
            torch.empty(rank, base.in_features, device=base.weight.device, dtype=base.weight.dtype)
        )
        self.lora_b = nn.Parameter(
            torch.zeros(base.out_features, rank, device=base.weight.device, dtype=base.weight.dtype)
        )
        nn.init.kaiming_uniform_(self.lora_a, a=math.sqrt(5))

    def forward(self, x):
        residual = (self.dropout(x).to(self.lora_a.dtype) @ self.lora_a.T) @ self.lora_b.T
        return self.base(x) + (self.scale * residual).to(x.dtype)


def inject_lora(model, rank=8, alpha=16, dropout=0.05):
    names = []
    for name, module in list(model.named_modules()):
        if (
            isinstance(module, nn.Linear)
            and name.startswith("transformer.blocks.")
            and name.endswith((".attn.layernorm_qkv.1", ".attn.out_proj"))
        ):
            parent_name, child_name = name.rsplit(".", 1)
            setattr(
                model.get_submodule(parent_name),
                child_name,
                LoRALinear(module, rank, alpha, dropout),
            )
            names.append(name)
    expected = 2 * len(model.transformer.blocks)
    if len(names) != expected:
        raise RuntimeError(
            f"ESM3 API changed: expected {expected} LoRA targets, found {len(names)}"
        )
    return names
