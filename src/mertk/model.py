"""Local-only ESM3 loading and adapter injection into its actual differentiable forward."""

import json
from pathlib import Path
from types import SimpleNamespace

import torch
from torch import nn
from torch.utils.checkpoint import checkpoint

from mertk.geometry import GeometricAdapter
from mertk.lora import inject_lora


def _unsupported(*args, **kwargs):
    raise RuntimeError(
        "Use this project's forward/design API; optional ESM3 decoders are not loaded"
    )


def make_esm3(tiny=False, device="cpu"):
    from accelerate import init_empty_weights
    from esm.models.esm3 import ESM3
    from esm.tokenization.sequence_tokenizer import EsmSequenceTokenizer

    # ESM3.forward only reads tokenizers.sequence. Functional annotation tokenizers
    # would unnecessarily download extra data. High-level encode/generate is not used.
    tokenizers = SimpleNamespace(sequence=EsmSequenceTokenizer())
    # Rotary inv_freq buffers are nonpersistent and absent from the state dict.
    # Keep buffers materialized while allocating the large parameters on meta.
    context = (
        init_empty_weights(include_buffers=False) if device == "meta" else torch.device(device)
    )
    with context:
        return ESM3(
            d_model=48 if tiny else 1536,
            n_heads=4 if tiny else 24,
            v_heads=4 if tiny else 256,
            n_layers=2 if tiny else 48,
            structure_encoder_fn=_unsupported,
            structure_decoder_fn=_unsupported,
            function_decoder_fn=_unsupported,
            tokenizers=tokenizers,
        )


def load_backbone(weights, device="cpu", tiny=False):
    if tiny:
        return make_esm3(tiny=True).to(device)
    weights = Path(weights)
    if not (weights / "manifest.json").exists():
        raise FileNotFoundError("Run scripts/download_weights.py to obtain verified local weights")
    model = make_esm3(device="meta")
    state = torch.load(
        weights / "data/weights/esm3_sm_open_v1.pth",
        weights_only=True,
        mmap=True,
        map_location="cpu",
    )
    model.load_state_dict(state, assign=True, strict=True)
    # Keep parameters fp32; autocast handles CUDA bf16 activations during training.
    return model.to(device=device, dtype=torch.float32)


def load_structure_encoder(weights, device="cpu"):
    from accelerate import init_empty_weights
    from esm.models.vqvae import StructureTokenEncoder

    with init_empty_weights(include_buffers=False):
        encoder = StructureTokenEncoder(
            d_model=1024, n_heads=1, v_heads=128, n_layers=2, d_out=128, n_codes=4096
        )
    state = torch.load(
        Path(weights) / "data/weights/esm3_structure_encoder_v0.pth",
        weights_only=True,
        mmap=True,
        map_location="cpu",
    )
    encoder.load_state_dict(state, assign=True, strict=True)
    return encoder.to(device=device, dtype=torch.float32).eval().requires_grad_(False)


def checkpoint_blocks(model):
    for block in model.transformer.blocks:
        original = block.forward

        def forward(*args, _original=original, **kwargs):
            if torch.is_grad_enabled():
                return checkpoint(_original, *args, use_reentrant=False, **kwargs)
            return _original(*args, **kwargs)

        block.forward = forward


class InterfaceModel(nn.Module):
    def __init__(
        self,
        backbone,
        use_lora=True,
        use_geometry=True,
        use_time=False,
        rank=8,
        alpha=16,
        dropout=0.05,
        hidden=128,
        gradient_checkpointing=True,
    ):
        super().__init__()
        self.backbone = backbone.requires_grad_(False)
        d_model = backbone.encoder.sequence_embed.embedding_dim
        self.lora_targets = inject_lora(backbone, rank, alpha, dropout) if use_lora else []
        self.adapter = GeometricAdapter(d_model, hidden) if use_geometry else None
        self.time_embed = (
            nn.Sequential(nn.Linear(3, hidden), nn.SiLU(), nn.Linear(hidden, d_model))
            if use_time
            else None
        )
        if self.time_embed is not None:
            nn.init.zeros_(self.time_embed[-1].weight)
            nn.init.zeros_(self.time_embed[-1].bias)
        if gradient_checkpointing:
            checkpoint_blocks(backbone)

    def forward(self, batch, sequence_tokens=None, time=None):
        tokens = batch["sequence_tokens"] if sequence_tokens is None else sequence_tokens
        condition = self.backbone.encoder.sequence_embed.weight.new_zeros(
            (*tokens.shape, self.backbone.encoder.sequence_embed.embedding_dim)
        )
        if self.adapter is not None:
            condition = condition + self.adapter(
                batch["backbone"],
                batch["neighbor_xyz"],
                batch["neighbor_features"],
                batch["neighbor_mask"],
            )
        if self.time_embed is not None:
            if time is None:
                raise ValueError("Time-conditioned model requires a noise level")
            t = time.reshape(-1, 1)
            encoding = torch.cat((t, torch.sin(torch.pi * t), torch.cos(torch.pi * t)), -1)
            condition = condition + self.time_embed(encoding).unsqueeze(1) * batch[
                "design_mask"
            ].unsqueeze(-1)
        # This hook applies once, before transformer blocks. It is removed before
        # checkpoint backward recomputation, which only replays the blocks.
        handle = self.backbone.encoder.register_forward_hook(
            lambda module, args, output: output + condition.to(output.dtype)
        )
        try:
            return self.backbone(
                sequence_tokens=tokens,
                structure_coords=batch["backbone"],
                structure_tokens=batch.get("structure_tokens"),
                chain_id=batch["chain_id"],
            ).sequence_logits
        finally:
            handle.remove()

    def adapter_state_dict(self):
        names = {name for name, parameter in self.named_parameters() if parameter.requires_grad}
        # Include geometry buffers too, while excluding all frozen base weights.
        return {
            name: value.detach().cpu()
            for name, value in self.state_dict().items()
            if name in names or name.startswith("adapter.")
        }

    def load_adapter_state_dict(self, state):
        expected = set(self.adapter_state_dict())
        if set(state) != expected:
            raise ValueError("Checkpoint adapter keys do not match the configured model")
        self.load_state_dict(state, strict=False)


def weight_identity(weights):
    return json.loads((Path(weights) / "manifest.json").read_text())
