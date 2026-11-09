# Baselines, ablations, and BindCraft controls

The [revised design brief](experiment_plan.md) prioritizes whether ESM3 plus its
adapter approaches ProteinMPNN on the same backbone pool **under shared AF3
assessment**, and whether the resulting candidates meet the MerTK protein-design
criteria. AF3 execution and structural scoring remain unimplemented. Recovery and
seed retention below are sequence diagnostics, not the final structural endpoint.

Run commands from the repository root with the project environment active.
Start with `configs/comparison_core.json`: frozen ESM3, LoRA, combined, and
ProteinMPNN. `configs/comparison_suite.json` retains all eight arms for follow-up
ablations. Both specify seeds 42/43/44, eight samples per
complex, temperature 1.0, and 16 reveal steps for the ESM3/atom-decoder sampler.
ProteinMPNN uses its own autoregressive sampler at the same temperature; these
are different generative models, not identical decoding algorithms.

## Held-out interface reconstruction

| Arm | Trainable parts | Information supplied |
| --- | --- | --- |
| `frozen_esm3` | None | Target/fixed binder sequences and both backbones/structure tokens |
| `lora` | LoRA | Same as frozen ESM3 |
| `geometry` | Atom adapter | ESM3 inputs plus target heavy-atom neighborhoods |
| `combined` | LoRA + atom adapter | Same as geometry |
| `combined_backbone_atoms` | LoRA + atom adapter | Same ESM3 inputs, only target backbone atoms retained in neighborhoods |
| `denoise` | LoRA + atom adapter + time embedding | Combined inputs with variable masking during training |
| `atom_decoder` | Atom adapter + small per-residue decoder | Same local target-atom features, no ESM3 or sequence-context model |
| `proteinmpnn` | None | Pretrained standard ProteinMPNN, both N/CA/C/O backbones and fixed sequences |

The backbone-only atom control filters the **same** up-to-32 target atom candidates;
it does not refill removed side chains with more distant backbone atoms. Target
residue identities remain available, including in backbone-atom feature labels.
Thus this contrast measures removal of side-chain atoms at a fixed neighborhood,
not removal of target sequence information. No adapter receives binder side chains.

```bash
python scripts/setup_external.py
# If examples predate the O-atom features, regenerate and recache them:
python scripts/prepare_data.py
bash scripts/cluster_data.sh
python scripts/cache_structure_tokens.py --device cuda

# Plan the core comparison, then execute/resume once the dataset is ready.
python scripts/run_comparisons.py --suite configs/comparison_core.json --output outputs/comparisons_core
python scripts/run_comparisons.py --suite configs/comparison_core.json --output outputs/comparisons_core --execute
python scripts/run_comparisons.py --suite configs/comparison_core.json --output outputs/comparisons_core --execute --resume

# Separate small implementation check; never interpret its metrics scientifically.
python scripts/run_comparisons.py --suite configs/comparison_core.json \
  --smoke-test --execute --output outputs/comparison_core_smoke

# Optional full suite remains available; the CLI default also uses this suite.
python scripts/run_comparisons.py --suite configs/comparison_suite.json --output outputs/comparisons_full
```

The core suite creates 12 arm/seed jobs, including six training runs. The full
suite creates 24 jobs, including 18 training runs. Training parameters
come from `configs/train.json` or `--base-config PATH`. The frozen ESM3 and
ProteinMPNN arms do not fine-tune. Checkpoints and results stay under the chosen
output directory. Use a new directory for a changed dataset or configuration.

`setup_external.py` downloads pinned upstream source and ProteinMPNN's
`vanilla_model_weights/v_48_020.pt` into ignored `external/`. Revisions are in
`configs/external_tools.json`; the loader rejects changed tracked source or weights.
ProteinMPNN runs in the project environment. BindCraft uses a separate environment.

Every arm receives identical split membership, per-complex design masks, and
per-example sampling seeds. All design positions have N/CA/C/O coordinates;
native design labels are explicitly removed before ProteinMPNN sampling. Target
and non-design binder residues must remain unchanged. Prepared file and mask
hashes are checked before pairing results. Resume validates the saved plan and
reuses completed results; interrupted training resumes from `last.pt`.

Each `comparison.json` stores sequences, input/mask identities, sample-average
native recovery, mean pairwise design-position Hamming diversity, and unique
sequence fraction. `summary.json` averages examples and seeds within each family,
then weights families equally. It includes paired recovery differences against
frozen ESM3 and, when present, ProteinMPNN, with 95% intervals from 2,000 family
bootstrap draws. These are sequence-metric contrasts; structural noninferiority
requires the separate planned AF3 endpoint. One family yields
no confidence interval; repeated assembly copies or seeds do not become
independent families. These intervals measure family variation after averaging
seeds, not a separate estimate of training-seed uncertainty.

The trainer's validation metrics use full-mask **single-pass argmax** recovery;
the comparison uses **sampled** recovery. Report these separately. Hyperparameters,
including each model's sampling temperature, should be tuned only on validation
families before freezing a final test protocol. The initial shared temperature
is a reproducible starting protocol, not a claim of optimal calibration.

The current 20-assembly pilot has no validation family. Planning works, but real
suite execution refuses it. Expand and audit the family-split data first.

## BindCraft original versus interface redesign

[BindCraft](https://github.com/martinpacesa/BindCraft) generates initial complexes.
Its accepted original binder is an **application control**, not a native PDB
sequence-recovery baseline. All redesign arms start from that same complex and
change only the same eligible interface positions. Keep the original control for
independent evaluation alongside every redesigned candidate.

On a Linux NVIDIA CUDA host with Conda installed:

```bash
python scripts/setup_external.py --tool bindcraft
# Uses the official installer in a separate BindCraft conda environment and
# downloads its AF2 parameters. Match CUDA to the allocated host.
bash scripts/install_bindcraft.sh --cuda 12.4 --pkg_manager conda

# Supply a selected, continuous extracellular-domain construct, not the whole receptor.
python scripts/prepare_bindcraft.py --target-pdb data/targets/mertk_domain.pdb \
  --target-chain A --name MerTK --lengths 65 120 --designs 10 \
  --output outputs/bindcraft_mertk
# Optional --hotspots A25,A30-33 uses the input PDB's author numbering.
bash outputs/bindcraft_mertk/run.sh

python scripts/import_bindcraft.py --input outputs/bindcraft_mertk/designs/Accepted \
  --binder-chain B --target-chain A --output data/application
python scripts/cache_structure_tokens.py --examples data/application/examples.jsonl --device cuda

# Uses existing trained checkpoints; never trains on these application examples.
python scripts/run_bindcraft_comparison.py --trained-suite outputs/comparisons_core/plan.json
python scripts/run_bindcraft_comparison.py --trained-suite outputs/comparisons_core/plan.json --execute
```

Create the job on the execution host: its generated paths are absolute. The job
copies its target PDB, upstream filters and advanced settings, records hashes and
residue mappings, and checks the maximum binder plus target fits the 512-token
budget. It refuses missing backbone atoms, ambiguous numbering gaps, nonexistent
hotspots, and nonempty job directories. The upstream environment is captured in
`external/BindCraft/environment-resolved.txt` and `pip-resolved.txt` after install;
the wrapper pins source, not every transitive installer dependency.

The importer retains `bindcraft_original.fasta` and source PDB provenance. It
marks all imported examples `application_holdout`; the application runner checks
exact sequence overlap with the training/validation/test manifest. Sequence-family
and domain-level TAM exclusion still depend on the training-set curation audit.
Original controls and all redesign sequences are exported together in
`outputs/bindcraft_comparison/candidates.jsonl`, with source complexes, chain IDs,
design positions, and `structure_evaluation: pending`.

Application summaries label similarity to the original as **seed retention**,
not native recovery. They report diversity without assigning a binding score.
Rescore original and redesigned candidates with the same independent complex
prediction/relaxation protocol before comparing interface confidence, backbone
consistency, clashes, and AXL/TYRO3 specificity. That structural scorer, glycan and
ligand-access checks, and experimental agonism assays are not implemented here.
Do not compare original BindCraft confidence scores directly with unrescored
redesigns or infer affinity from sequence retention.

## Verified setup

The eight-arm synthetic suite has completed with tiny actual ESM3 modules and
the downloaded standard [ProteinMPNN](https://github.com/dauparas/ProteinMPNN)
weights. ProteinMPNN also sampled a real prepared pilot interface while preserving
fixed positions. These are integration checks, not benchmark results. BindCraft
source and example job preparation were checked; CUDA environment installation
and MerTK generation have not run.
