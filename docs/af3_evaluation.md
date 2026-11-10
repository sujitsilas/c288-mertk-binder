# Shared AF3 structural evaluation

Both ProteinMPNN and ESM3 variants redesign the same interfaces. This workflow
assesses their sequences with the same AF3 inputs and prediction budget. It
measures **predicted structural consistency**. Every output retains
`biological_evaluation: not_assessed`; a structural pass does not establish affinity,
specificity, ligand compatibility or efferocytosis.

## Train and export matched candidates

```bash
python scripts/run_comparisons.py                    # Core suite, plan only
python scripts/run_comparisons.py --execute          # Requires adequate train/val/test families
python scripts/run_comparisons.py --execute --resume
```

Comparisons use `best.pt`, selected by validation family-averaged reconstruction
recovery. Resume uses `last.pt`. The supervised objective remains masked sequence
reconstruction; AF3 is a downstream evaluator and is not a training loss.

Real runs export `outputs/comparisons_core/candidates.jsonl`, including an original
sequence control per complex. The BindCraft application runner produces the same
schema in `outputs/bindcraft_comparison/candidates.jsonl`. Candidate records include
family/split, full sequences, design positions, reference NPZ/mask hashes and model
provenance. Each non-control arm must contain identical complex/seed/sample keys.
Synthetic training output is refused as structural benchmark input.

## Prepare folding jobs

`configs/af3_evaluation.json` starts in **draft** state. It specifies the required
four arms, AF3 seed 101, five diffusion samples per candidate, ten recycles, and the
structural thresholds described in the experiment plan. AF3's diffusion sampling
is independent of the optional ESM sequence-denoising training objective.

The primary protocol uses a cached target MSA, identical across arms and binder
variants. Prepare it using the AF3 data pipeline or your existing MSA workflow and
record its database version. Create a JSON mapping the SHA-256 of each complete
target sequence to its A3M file; paths resolve relative to that mapping file:

```json
{
  "FULL_SHA256_OF_TARGET_SEQUENCE": "mertk_target.a3m"
}
```

The preparer prints the needed sequence hash if one is missing. The first A3M
record must exactly match the target; other aligned rows must have matching length
after removal of lowercase insertions. Both chains use empty structural template
lists. The binder has no MSA. The designed complex is retained only as a scoring
reference and is never passed to AF3 as a template.

```bash
python scripts/prepare_af3.py \
  --candidates outputs/comparisons_core/candidates.jsonl \
  --target-msas data/target_msas/index.json --output outputs/af3_core
```

For an explicitly MSA-free exploratory protocol, use `--single-sequence` instead
of `--target-msas`. This choice is recorded; it is not substituted silently for
missing target MSAs. A locked protocol cannot be overridden with that flag.

To bound an integration pilot, add `--max-complexes 6 --samples-per-design-seed 4`.
Selection is deterministic and independent of model scores: complex IDs are
ordered by their hash, and the first requested sample indices are retained for
each method/seed. Original controls are preserved. This limits candidates, not
independent families; inspect the selected cohort before committing resources.

The output bundle contains `inputs/`, copied `references/`, and `manifest.json`.
It is portable as a directory. Input hashes, reference hashes, protocol and cohort
are immutable; changed settings require a new bundle. The printed prediction count
includes every candidate, original control, AF3 seed and diffusion sample.

## Execute in a separate AF3 environment

```bash
python scripts/setup_external.py --tool alphafold3
```

This downloads source only. Follow the pinned upstream
[installation instructions](https://github.com/google-deepmind/alphafold3/blob/28c2d290559fe07387a12f471ad49f0d36206327/docs/installation.md)
to prepare the AF3 runtime and model parameters on the inference host. The project
environment does not install AF3/JAX or obtain AF3 parameters automatically.

From the project environment, point the launcher at that separate interpreter:

```bash
# Review the command without launching inference.
python scripts/run_af3.py --manifest outputs/af3_core/manifest.json \
  --af3-python /path/to/af3-env/bin/python --model-dir /path/to/af3-parameters

# Run all pending candidates in one AF3 process; add --resume after interruption.
python scripts/run_af3.py --manifest outputs/af3_core/manifest.json \
  --af3-python /path/to/af3-env/bin/python --model-dir /path/to/af3-parameters --execute
```

The launcher requires the pinned, unmodified source checkout, hashes parameter
files, and records Python/package versions. Resume requires the same engine and
plan. It passes precomputed MSAs with the AF3 data pipeline disabled and preserves
the common seed/sample/recycle settings. Incomplete jobs can resume; completed jobs
with matching processed inputs are skipped. Wall time is recorded for the whole
pending batch, including startup/compilation; it is not a per-arm AF3 timing.

## Score and inspect results

```bash
python scripts/score_af3.py --manifest outputs/af3_core/manifest.json \
  --predictions outputs/af3_core/predictions --output outputs/af3_core_scores
```

The scorer checks the processed AF3 sequences/MSAs/templates, every planned
seed/sample, chain identities, residue mappings and finite confidence values.
It selects the highest-ranked prediction from the full-precision AF3 ranking CSV,
using upstream seed/sample order to resolve exact ties. It never selects the
prediction with the smallest reference RMSD. Summary confidence files round their
scores, so selecting directly from those files could choose a different sample.

For target A and binder B it reports:

- Pair ipTM and the AF3 severe-clash flag (numeric 0/1 or boolean).
- Binder Cα RMSD after fitting the binder, and binder Cα RMSD after fitting the
  target. The second measurement detects an incorrect binding pose despite a
  correct binder fold.
- Target fit RMSD and reference-coordinate coverage. Alignment uses all finite
  reference Cα positions in the prepared construct; collinear fits are rejected.
- Retention of reference interchain Cα contacts at 8 Å, directional whole-chain
  mean PAE, and mean binder atom pLDDT. The Cα contact diagnostic differs from the
  CA-to-target-heavy-atom mask used for sequence design.

The draft structural gate requires pair ipTM ≥0.8, binder fold RMSD ≤2 Å, target-
aligned binder RMSD ≤3 Å, at least 90% reference Cα coverage in each chain, and no
severe AF3 clash. These are starting engineering cutoffs, not validated binding-hit
rules. Inspect detailed local clashes separately; the global clash flag can miss them.

`candidates.jsonl` retains every planned candidate and selected-file hashes.
Missing, malformed or mismatched predictions are explicitly marked. If any planned
output is incomplete, `summary.json` reports counts and valid-output rates but
withholds paired estimates for the cohort. Failed structural gates on otherwise
valid predictions count as failures and remain in the denominator.

Complete cohorts receive family-averaged pass rates and paired differences versus
ProteinMPNN, LoRA and frozen ESM3, with 2,000 family bootstrap draws. Original
controls are reported separately. Seeds, samples and repeated complexes are not
independent families. The optional noninferiority flag requires a locked protocol,
verified launcher provenance, a complete test cohort and the configured minimum
number of families (default ten); it uses the lower bootstrap bound and a 0.05
absolute margin. That minimum is a guard, not a power calculation. Bootstrap
intervals describe sampled-family variation after averaging seeds and can be
unreliable with few families. Draft, validation and application runs never assert
noninferiority. Calibrate thresholds on validation, justify the margin and budget,
then copy the protocol with `protocol_state: locked` before final test execution.

## Validation status

Unit/integration fixtures cover input matching, cached-MSA checks, fixed-position
preservation, rigid-motion invariance, fold versus pose, AF3 ranking precision,
missing/invalid predictions, and inference guards. The source revision and actual
upstream JSON serialization were inspected. No AF3 model inference or real
structural benchmark has been run in this repository yet.
