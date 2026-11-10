# MerTK interface sequence design

The protein-design goal is a MerTK-selective extracellular binder that preserves
natural-ligand access, as a building block for future tests of macrophage
efferocytosis. The computational question is whether **ESM3 with LoRA and a
target-atom adapter can approach ProteinMPNN for interface sequence redesign,
with both evaluated through the same AF3 complex-prediction protocol**.

See the [revised project direction](Report/revised_project_direction.md), responding
to feedback on the [original proposal](Report/final_project_proposal.pdf), and the
[design criteria and experiment plan](docs/experiment_plan.md). Native sequence
reconstruction trains the model; it is a diagnostic rather than the final success
criterion. MerTK, AXL, and TYRO3 remain reserved for the application stage.

Sequence training, sampling, matched AF3 job preparation/launching, and structural
pass-rate scoring are implemented. AF3 integration has been tested with output
fixtures; real AF3 inference and biological compatibility checks have not run.
Improved binding, receptor agonism, and efferocytosis have not been established.

## Environment and weights

Run commands from the repository root:

```bash
uv sync --frozen --extra train --extra dev
source .venv/bin/activate
export MPLCONFIGDIR="$PWD/.cache/matplotlib"
python scripts/download_weights.py --dry-run
python scripts/download_weights.py
```

Python 3.12 and `esm==3.4.2` are pinned. `configs/weights.json` pins the upstream
checkpoint revision. The downloader retrieves the ESM3-open 1.4B model and its
backbone structure encoder, verifies sizes and SHA-256 checksums, and writes a
local manifest. Downloads resume through Hugging Face's cache. Use `HF_TOKEN` or
`hf auth login` if your environment requires authentication; never put tokens in git.
Weights, datasets, environments, and training outputs are ignored by git.

The upstream model moved to [Biohub](https://huggingface.co/biohub/esm3-sm-open-v1).
The official [`esm` implementation](https://github.com/Biohub/esm) provides the
model API. No AlphaFold 3 or RFD4 weights are required for adapter training.
The AF3 assessment needs a separate installed runtime and model parameters;
see [the structural evaluation guide](docs/af3_evaluation.md).

The selected files total **2.87 GB**. On 2026-11-04 these were downloaded and
verified locally. Downloads are not in the GitHub repository; run the script on
each training machine or transfer your ignored `weights/` directory.

## Check the implementation

```bash
pytest -q
python scripts/train.py --smoke-test
python scripts/train.py --smoke-test --objective denoise --output outputs/smoke_denoise
# Optional: uses the actual 1.4B weights, several GB RAM, and no network.
python scripts/check_pretrained.py --backward
```

Smoke tests use tiny, randomly initialized **actual ESM3 modules** with synthetic
data. They check implementation, not model quality. Choose a new `--output` for
repeated smoke runs. [Validation notes](docs/validation.md) distinguish these from
the completed full-checkpoint integration check.

## Prepare training data

Install [MMseqs2](https://github.com/soedinglab/MMseqs2#installation) separately
(`brew install mmseqs2` on macOS, or use the official Linux/Conda installation).
Python dependencies alone do not install its executable.

```bash
# Bounded download; increase after validating the pilot.
python scripts/fetch_pdb.py --limit 2000
python scripts/prepare_data.py --max-length 512
bash scripts/cluster_data.sh
python scripts/cache_structure_tokens.py --device cuda
```

Use `--device cpu` for structure token caching on a laptop. PDB discovery uses the
[RCSB Search API](https://search.rcsb.org/), followed by expanded biological
assembly mmCIF downloads. It requires protein-only assemblies, at least two protein
chains, X-ray/cryo-EM, and resolution ≤3 Å. The number requested is the number of
assemblies **before** interface/length/homomer filters. `--assemblies 1BRS-1` can
fetch an explicit assembly for debugging, bypassing search quality filters.

Preprocessing writes full polymer sequences (including missing-coordinate
positions), backbone coordinates, target-only atom neighborhoods, fixed design
masks, and provenance. Chains longer than the joint 512-token budget are skipped;
there is no silent sequence truncation. The default contact rule is binder CA to
target heavy atom ≤8 Å, with at least five interface residues and complete
N/CA/C/O coordinates at designed positions. This coarse rule
avoids defining masks from the native binder side chains. Neighborhoods contain
up to 32 target atoms within 12 Å. Only N/CA/C enter the frozen structure tokenizer.

MMseqs2 clusters the **binder and target** sequences plus UniProt TAM full-sequence
and domain seeds at 30% identity / 80% bidirectional coverage. Connected components
join all chain families and all assemblies of a PDB entry before deterministic
80/10/10 train/validation/test assignment. Components touching a TAM seed are
reserved as `application_holdout`. Same-family binder/target pairs are excluded
because the visible target can reveal the answer. Read the limitations of this
sequence-based family approximation in [the experiment plan](docs/experiment_plan.md).

The initial local pilot contains 20 downloaded assemblies and 16 directed
interfaces from two PDB entries. Its split has 8 train and 8 test examples and
**no validation families**. Training deliberately refuses this undersized pilot.
Expand the data; do not move related chains between folds to fill a missing split.
The data and processed features are local and ignored by git.

## Train and compare

Start with the [four-arm core comparison](configs/comparison_core.json):
ProteinMPNN, frozen ESM3, LoRA, and LoRA plus the geometric adapter. It uses three
seeds and common design masks. Current summaries report sampled recovery/diversity
and paired family-level differences versus ProteinMPNN and frozen ESM3. These
sequence statistics remain separate from the [AF3 structural endpoint](docs/af3_evaluation.md).

The [comparison guide](docs/comparisons.md) also retains the full eight-arm suite
for adapter-only, atom-only decoder, backbone-only atom, and denoising follow-ups.

```bash
python scripts/setup_external.py  # Pinned ProteinMPNN weights and BindCraft source
python scripts/run_comparisons.py --suite configs/comparison_core.json \
  --output outputs/comparisons_core  # Plan 12 jobs; six involve training
python scripts/run_comparisons.py --suite configs/comparison_core.json \
  --output outputs/comparisons_core --execute
# Complete implementation check with tiny ESM3 modules and actual ProteinMPNN:
python scripts/run_comparisons.py --suite configs/comparison_core.json \
  --smoke-test --execute --output outputs/comparison_core_smoke
```

BindCraft has a separate CUDA environment installer, job preparer, accepted-complex
importer, and application comparison runner. Original BindCraft sequences are
retained as controls alongside redesigns; see the guide for commands. Complex
structure evaluation uses the shared AF3 workflow, and application seed retention
is not treated as native recovery or improved binding.

Run in an allocated CUDA GPU session on Hoffman2 or a lab workstation. The current
trainer uses one GPU, bf16 autocast, fp32 parameters, gradient checkpointing,
one complex per microbatch, and four accumulation steps. Start with a 48 GB RTX
6000 Ada; actual CUDA memory use and throughput still need measurement.
Two GPUs can run separate ablations. There is no distributed trainer or automatic
cluster job submission.

```bash
python scripts/train.py --ablation combined --output outputs/combined
python scripts/train.py --ablation lora --output outputs/lora
python scripts/train.py --ablation geometry --output outputs/geometry
# Optional follow-up after the core comparison:
python scripts/train.py --objective denoise --output outputs/denoise

# Resume an interrupted run with the same config, objective, and ablation.
python scripts/train.py --output outputs/combined --resume outputs/combined/last.pt

# Same full-interface masks for every evaluation.
python scripts/evaluate.py --split test --output outputs/frozen_test.json
python scripts/evaluate.py --checkpoint outputs/combined/best.pt \
  --split test --output outputs/combined_test.json
```

`configs/train.json` controls learning rate, LoRA rank, accumulation, and run
length. The default combined setup has **3,758,321 trainable parameters** (~0.27%
of the complete model). The backbone remains frozen; LoRA wraps the fused QKV
and output projections in all 48 attention blocks. Checkpoints contain adapters,
optimizer/RNG state, base-weight identity, configuration, git revision, and split
manifest hash. They do not duplicate ESM3's pretrained weights. Metrics report
masked-token loss and 20-amino-acid interface recovery.

Training samples families uniformly, then examples within each family. Validation
also reports family-averaged recovery. `best.pt` is selected by that validation
metric (earliest step wins ties); comparisons use it. `last.pt` retains optimizer
and RNG states for exact resume. Both are adapter-only checkpoints. Use
`--train-family-fraction` on the comparison runner for matched nested data subsets.
These checkpoints optimize a reconstruction diagnostic, not AF3 confidence or affinity.

The geometric adapter expresses target atoms in each binder residue's N/CA/C
frame, embeds element/residue/backbone features and radial distances, pools them
with learned attention, and adds a zero-initialized projection before the ESM3
transformer. Both target sequence and backbone remain visible in **all ESM3** ablations.
Frozen structure tokens are cached separately for each chain.

`--objective denoise` adds a noise-level embedding and trains reconstruction under
variable masking rates. The sampler progressively fills masked residues with an
absorbing-mask reverse process. This is an optional **sequence denoising**
experiment; it does not implement coordinate diffusion, RFD4 latent flow matching,
or the AF3 Pairformer. Its cross-entropy objective is not a diffusion likelihood
bound. Those larger extensions should follow a successful reconstruction baseline.

## Sample fixed-backbone sequences

```bash
python scripts/design.py --checkpoint outputs/combined/best.pt \
  --example data/processed/ASSEMBLY_BINDER_TARGET.npz \
  --samples 8 --steps 16 --output outputs/designs.fasta
```

Replace the example path with a prepared interface whose structure tokens have
been cached. The target, backbone, and non-design residues remain fixed. FASTA and
a provenance JSON record are written. Unique-sequence count is reported; binding,
fold consistency, glycan clashes, and ligand compatibility are not scored yet.

See [the proposal review and experiment plan](docs/experiment_plan.md) for the
AF3/RFD4 connection, benchmark controls, and remaining application work.
