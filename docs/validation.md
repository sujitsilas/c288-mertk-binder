# Validation record — 2026-11-04

Validated locally on macOS/CPU using Python 3.12.11, ESM 3.4.2, and PyTorch 2.11.0.
The exact Python environment is recorded in `uv.lock`.

## Completed checks

- Downloaded the pinned ESM3-open base model and structure encoder (2.87 GB),
  checked upstream file sizes and SHA-256 values, and wrote `weights/esm3/manifest.json`.
- Retrieved 20 RCSB biological assemblies and the three UniProt TAM records.
  Preprocessing yielded 16 directed interfaces from two PDB entries.
- Encoded all 16 examples with the actual pretrained ESM3 structure encoder.
- Ran MMseqs2 connected-component clustering and split assignment. The pilot has
  8 training and 8 test examples, with no validation examples. This is a pipeline
  fixture, not a usable training/evaluation benchmark.
- Ran masked-reconstruction and time-conditioned denoising training smoke checks
  with a tiny random ESM3, including validation and adapter-only checkpoint writes.
- Ran the **actual 1.4B pretrained model** forward and backward on a synthetic
  25-token tensor fixture, checking nonzero LoRA/geometry gradients and no frozen
  backbone gradients. Output shape was `[1, 25, 64]`, with 96 LoRA insertion sites,
  3,758,321 trainable parameters, and 1,405,494,069 total parameters.
- Unit tests cover rigid-motion invariance, missing-frame/empty-neighbor behavior,
  binder-side-chain input isolation, family/entry/TAM split closure, cross-fold
  sequence leakage, adapter initialization equivalence, ESM3 gradients, local
  checkpoint loading, checkpoint round trips, deterministic optimizer resume,
  and fixed-residue sequence sampling.

```bash
pytest -q
ruff check src scripts tests
ruff format --check src scripts tests
```

The full-model integration report is local at `outputs/pretrained_check.json`.
All training-like checks used synthetic data. No substantive fine-tuning run,
held-out performance improvement, CUDA memory claim, or MerTK binder result has
been established. No cluster jobs were submitted.

## Integration issues caught and fixed

ESM3's rotary buffers are absent from its checkpoint state dictionary. Loading
parameters on the meta device must leave those nonpersistent buffers materialized;
the loader uses `accelerate.init_empty_weights(include_buffers=False)`.

LoRA target names explicitly match `.attn.*` so they do not accidentally match
the separate `.geom_attn.out_proj` projection.

The tested MMseqs2 build requires `--single-step-clustering 1` when using connected
component mode with the current default linclust implementation. Its version was
`2a3dafc46f9a1bc79af96b4c11a24d63a7b0b4d4`; the pipeline records the installed
version for each future run. MMseqs2 was run from a temporary portable installation;
install it on the training machine using the README instructions.

## Reoriented comparison checks — 2026-11-09

- Added the four-arm core configuration: frozen ESM3, LoRA, combined, and
  ProteinMPNN. Planning produced 12 arm/seed jobs, including six training runs.
- Completed the four-arm synthetic smoke suite with tiny ESM3 modules and the
  actual downloaded ProteinMPNN weights. Its report is local at
  `outputs/comparison_core_smoke_20261109/summary.json`.
- Sequence summaries now include paired family-level differences versus
  ProteinMPNN as well as frozen ESM3. Tests cover both native recovery and
  application seed-retention labels, including mismatched-cohort rejection.
- `pytest -q`: 36 passed, 1 skipped. Ruff checks and formatting checks passed.

At the reorientation commit, the structural pass-rate endpoint was a proposed
protocol. The subsequent implementation is recorded below; the current pilot still
lacks independent validation families. No scientific performance comparison or
new binder-function result follows from these implementation checks.

## Training and structural evaluator implementation — 2026-11-10

- The four-arm suite is now the default. Family-balanced sampling, nested training
  family subsets, separate data/noise RNG streams and validation-family checkpoint
  selection are implemented. Both ordinary and accelerator resume checks passed.
- The updated four-arm training/sampling smoke suite completed and resumed from
  its existing results. Reports are in `outputs/core_training_v2_smoke/` locally.
  Smoke results are blocked from becoming structural benchmark candidates.
- Downloaded and inspected AF3 source pinned at
  `28c2d290559fe07387a12f471ad49f0d36206327`. The actual serializer uses a numeric
  clash flag and rounds summary scores; the scorer accepts the flag and selects
  samples using the full-precision ranking CSV.
- Fixture tests cover candidate export, immutable job inputs, target MSA matching,
  fixed-position preservation, explicit chain/residue mapping, rigid transforms,
  binder fold versus target-relative pose, ranking selection, missing/invalid
  outputs, paired family comparisons and noninferiority guards.
- `pytest -q`: 53 passed, 1 skipped. The accelerator resume test passed separately
  outside the restricted sandbox. Ruff checks and formatting checks passed.

AF3 model weights/environment have not been installed here and no real AF3 inference
has run. Biological compatibility and functional assays remain unassessed. The
structural thresholds and margin remain a draft protocol pending validation.
