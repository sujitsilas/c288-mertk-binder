#!/usr/bin/env python3
"""Exercise a real downloaded ESM3 checkpoint on a small synthetic tensor fixture."""

import argparse
import json
from pathlib import Path

import torch

from mertk.dataset import make_batch, synthetic_example
from mertk.model import InterfaceModel, load_backbone
from mertk.objectives import corrupt, reconstruction_loss


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--weights", default="weights/esm3")
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--backward", action="store_true")
    parser.add_argument("--output", type=Path, default=Path("outputs/pretrained_check.json"))
    args = parser.parse_args()
    torch.set_num_threads(4)
    torch.manual_seed(42)
    base = load_backbone(args.weights, args.device)
    model = InterfaceModel(base, dropout=0).to(args.device)
    batch = make_batch(synthetic_example(), base.tokenizers.sequence, args.device)
    tokens, selected = corrupt(batch["sequence_tokens"], batch["design_mask"], 32, 1.0)
    model.eval()
    with torch.set_grad_enabled(args.backward):
        logits = model(batch, sequence_tokens=tokens)
        loss = reconstruction_loss(logits, batch["sequence_tokens"], selected)
        if not torch.isfinite(loss):
            raise FloatingPointError("Non-finite pretrained model loss")
        if args.backward:
            loss.backward()
            assert model.adapter.out.weight.grad.abs().sum() > 0
            assert any(
                p.grad is not None and p.grad.abs().sum() > 0
                for name, p in model.named_parameters()
                if "lora_b" in name
            )
            assert all(p.grad is None for p in model.parameters() if not p.requires_grad)
    result = {
        "pretrained": True,
        "synthetic_input": True,
        "scientific_benchmark": False,
        "device": args.device,
        "backward_verified": args.backward,
        "loss": loss.item(),
        "logits_shape": list(logits.shape),
        "lora_targets": len(model.lora_targets),
        "trainable_parameters": sum(p.numel() for p in model.parameters() if p.requires_grad),
        "total_parameters": sum(p.numel() for p in model.parameters()),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
