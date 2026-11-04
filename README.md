# MerTK interface sequence design

Implementation of [the project proposal](Report/final_project_proposal.pdf):
fine-tune ESM3-open with LoRA and a target-atom geometric adapter, then test whether
it improves reconstruction of interface sequences on held-out protein families.
MerTK, AXL, and TYRO3 are reserved for the application stage.

This is a research starting point. Improved binding, receptor agonism, and
efferocytosis are hypotheses, not established outcomes of sequence reconstruction.

## Environment and weights

```bash
uv sync --extra train --extra dev
uv run python scripts/download_weights.py --dry-run
uv run python scripts/download_weights.py
```

Python 3.12 and `esm==3.4.2` are pinned. `configs/weights.json` pins the upstream
checkpoint revision. The downloader retrieves the ESM3-open 1.4B model and its
backbone structure encoder, verifies sizes and SHA-256 checksums, and writes a
local manifest. Downloads resume through Hugging Face's cache. Use `HF_TOKEN` or
`hf auth login` if your environment requires authentication; never put tokens in git.
Weights, datasets, environments, and training outputs are ignored by git.

The upstream model moved to [Biohub](https://huggingface.co/biohub/esm3-sm-open-v1).
The official [`esm` implementation](https://github.com/Biohub/esm) provides the
model API. No AlphaFold 3 or RFD4 weights are required for this adapter experiment.

Pipeline commands, architecture, and validation details are added alongside their
implementation in subsequent commits.
