"""Masked reconstruction and optional time-conditioned absorbing-mask denoising."""

import torch
from torch.nn import functional as F


def corrupt(tokens, design_mask, mask_token_id, probability=0.5, generator=None):
    if not 0 < probability <= 1:
        raise ValueError("Mask probability must lie in (0, 1]")
    if not design_mask.any(dim=1).all():
        raise ValueError("Every complex must contain designable positions")
    # CPU generator makes masks identical across CPU/GPU evaluation and ablations.
    random = torch.rand(tokens.shape, generator=generator).to(tokens.device)
    selected = (random < probability) & design_mask
    for i in range(len(tokens)):
        if not selected[i].any():
            candidates = design_mask[i].nonzero().flatten()
            j = int(torch.randint(len(candidates), (), generator=generator))
            selected[i, candidates[j]] = True
    corrupted = tokens.clone()
    corrupted[selected] = mask_token_id
    return corrupted, selected


def reconstruction_loss(logits, labels, selected):
    if not selected.any():
        raise ValueError("Cannot compute reconstruction loss on an empty mask")
    # Use the full ESM vocabulary for loss; restrict design generation to the 20 AAs.
    return F.cross_entropy(logits[selected].float(), labels[selected])


@torch.no_grad()
def sample_sequence(model, batch, tokenizer, steps=16, temperature=1.0, seed=42):
    """Absorbing-mask reverse sampler: reveal masks with probability 1/remaining_steps.

    Coordinates and non-design tokens stay fixed. This is a discrete sequence
    denoising experiment, not an AF3/RFD4 coordinate or latent diffusion model.
    """
    if steps < 1 or temperature <= 0:
        raise ValueError("steps and temperature must be positive")
    model.eval()
    tokens = batch["sequence_tokens"].clone()
    design_mask = batch["design_mask"]
    tokens[design_mask] = tokenizer.mask_token_id
    generator = torch.Generator().manual_seed(seed)
    from mertk.data import AA

    aa_ids = torch.tensor([tokenizer.convert_tokens_to_ids(a) for a in AA], device=tokens.device)
    for remaining in range(steps, 0, -1):
        time = torch.full((len(tokens),), remaining / steps, device=tokens.device)
        logits = model(batch, sequence_tokens=tokens, time=time)
        pending = design_mask & (tokens == tokenizer.mask_token_id)
        if not pending.any():
            break
        probabilities = (logits[pending][:, aa_ids].float() / temperature).softmax(-1).cpu()
        choices = torch.multinomial(probabilities, 1, generator=generator).squeeze(-1)
        reveal = torch.rand(len(choices), generator=generator) < 1 / remaining
        pending_indices = pending.nonzero()
        chosen_indices = pending_indices[reveal.to(tokens.device)]
        tokens[chosen_indices[:, 0], chosen_indices[:, 1]] = aa_ids[
            choices[reveal].to(tokens.device)
        ]
    if (tokens[design_mask] == tokenizer.mask_token_id).any():
        raise RuntimeError("Sampler failed to fill all design positions")
    return tokens
