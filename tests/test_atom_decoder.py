import torch

from mertk.atom_decoder import AtomDecoder
from mertk.dataset import make_batch, synthetic_example
from mertk.objectives import corrupt, reconstruction_loss


def test_atom_decoder_trains_and_roundtrips_without_an_esm_backbone():
    model = AtomDecoder(hidden=32)
    batch = make_batch(synthetic_example(), model.tokenizer)
    tokens, selected = corrupt(batch["sequence_tokens"], batch["design_mask"], 32, 1.0)
    optimizer = torch.optim.AdamW(model.parameters(), lr=0.001)
    loss = reconstruction_loss(
        model(batch, sequence_tokens=tokens), batch["sequence_tokens"], selected
    )
    loss.backward()
    assert model.adapter.out.weight.grad.abs().sum() > 0
    optimizer.step()
    restored = AtomDecoder(hidden=32)
    restored.load_adapter_state_dict(model.adapter_state_dict())
    torch.testing.assert_close(restored(batch), model(batch))
    assert not any(name.startswith("backbone") for name, _ in model.named_parameters())
