import torch

from mertk.dataset import make_batch, synthetic_example
from mertk.model import InterfaceModel, load_backbone, make_esm3
from mertk.objectives import corrupt, reconstruction_loss, sample_sequence


def setup_model(**kwargs):
    torch.manual_seed(7)
    base = make_esm3(tiny=True)
    model = InterfaceModel(base, dropout=0, **kwargs)
    batch = make_batch(synthetic_example(), base.tokenizers.sequence)
    return model, batch


def test_zero_initialized_adapters_preserve_actual_esm3_forward():
    base = make_esm3(tiny=True).eval()
    batch = make_batch(synthetic_example(), base.tokenizers.sequence)
    with torch.no_grad():
        expected = base(
            sequence_tokens=batch["sequence_tokens"],
            structure_coords=batch["backbone"],
            chain_id=batch["chain_id"],
        ).sequence_logits
    model = InterfaceModel(base, dropout=0).eval()
    with torch.no_grad():
        torch.testing.assert_close(model(batch), expected)
    assert len(model.lora_targets) == 4


def test_local_loader_materializes_rotary_buffers(tmp_path, monkeypatch):
    model = make_esm3(tiny=True)
    (tmp_path / "data/weights").mkdir(parents=True)
    (tmp_path / "manifest.json").write_text("{}")
    torch.save(model.state_dict(), tmp_path / "data/weights/esm3_sm_open_v1.pth")
    monkeypatch.setattr("mertk.model.make_esm3", lambda device: make_esm3(tiny=True, device=device))
    loaded = load_backbone(tmp_path)
    assert not any(buffer.is_meta for buffer in loaded.buffers())
    for name, parameter in loaded.named_parameters():
        torch.testing.assert_close(parameter, model.get_parameter(name))


def test_real_esm3_backward_updates_only_adapters_with_checkpointing():
    model, batch = setup_model(gradient_checkpointing=True)
    original = {
        name: p.detach().clone() for name, p in model.named_parameters() if not p.requires_grad
    }
    trainable = [p for p in model.parameters() if p.requires_grad]
    optimizer = torch.optim.AdamW(trainable, lr=0.01)
    tokens, selected = corrupt(batch["sequence_tokens"], batch["design_mask"], 32, 1.0)
    for _ in range(2):
        optimizer.zero_grad(set_to_none=True)
        loss = reconstruction_loss(
            model(batch, sequence_tokens=tokens), batch["sequence_tokens"], selected
        )
        loss.backward()
        optimizer.step()
    assert model.adapter.atom_mlp[0].weight.grad.abs().sum() > 0
    assert any(
        p.grad is not None and p.grad.abs().sum() > 0
        for name, p in model.named_parameters()
        if "lora_a" in name
    )
    for name, p in model.named_parameters():
        if not p.requires_grad:
            assert p.grad is None
            torch.testing.assert_close(p, original[name], atol=0, rtol=0)


def test_checkpoint_roundtrip_and_no_frozen_weights(tmp_path):
    model, batch = setup_model(use_time=True)
    with torch.no_grad():
        model.adapter.out.weight.normal_()
    path = tmp_path / "adapter.pt"
    torch.save(model.adapter_state_dict(), path)
    loaded = torch.load(path, weights_only=True)
    assert not any("base.weight" in key or "sequence_embed" in key for key in loaded)
    restored, _ = setup_model(use_time=True)
    restored.load_adapter_state_dict(loaded)
    model.eval()
    restored.eval()
    time = torch.ones(1)
    torch.testing.assert_close(model(batch, time=time), restored(batch, time=time))


def test_corruption_and_sampling_keep_target_and_special_tokens_fixed():
    model, batch = setup_model(use_time=True)
    tokenizer = model.backbone.tokenizers.sequence
    original = batch["sequence_tokens"].clone()
    fixed = ~batch["design_mask"]
    corrupted, selected = corrupt(
        original,
        batch["design_mask"],
        32,
        probability=0.001,
        generator=torch.Generator().manual_seed(1),
    )
    assert selected.any()
    assert torch.equal(corrupted[fixed], original[fixed])
    generated = sample_sequence(model, batch, tokenizer, steps=3, seed=7)
    assert torch.equal(generated[fixed], original[fixed])
    assert ((generated[batch["design_mask"]] >= 4) & (generated[batch["design_mask"]] < 24)).all()
    assert torch.equal(original, batch["sequence_tokens"])
    assert torch.equal(generated, sample_sequence(model, batch, tokenizer, steps=3, seed=7))
