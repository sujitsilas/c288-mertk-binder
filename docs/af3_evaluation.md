# Shared AF3 assessment of externally designed complexes

AF3 evaluates supplied candidate sequences after design. The receptor starting reference is experimental 37NT; AF3 predictions are a consistency check and do not replace it. This implementation handles **one target protein chain A and one binder protein chain B**, with no glycans, membranes or ligand chains in the scoring protocol. Separate biological review and specialized follow-up complexes are required.

## Prepare the candidate cohort

Import completed designs and merge their manifests as described in [comparisons](comparisons.md). Every arm must cover the same campaign/epitope strata with matching target sequences and experimental-reference checksums. Each candidate has its own design-complex reference. Different generators may produce different binder backbones, so there is no requirement that their complex IDs or sequences be paired.

The reference NPZ contains the imported sequence and backbone coordinates. Hash checks reject changes after import; sequence checks prevent scoring a reference under a different candidate sequence. Importing a structure records it as pending, not as a successful binder. Experimental-reference identity does not certify that the generator respected that structure: audit target sequence, coordinate distortion, crop mapping and glycans separately.

## Protocol and job preparation

The draft [protocol](../configs/af3_evaluation.json) uses one model seed, five predictions per seed, ten recycles, and the same target MSA policy for every arm. These are starting settings; confirm a feasible runtime budget and freeze the protocol before comparing nominated designs.

The default requires cached target A3M files. Supply a JSON object mapping the SHA-256 hash of each target sequence to its A3M path, resolved relative to that JSON file. The A3M query must exactly match the target construct. The binder uses a single sequence. Neither protein has structural templates in this assessment, to avoid imposing the proposed binder pose.

```bash
python scripts/prepare_af3.py \
  --candidates data/candidates/pilot.jsonl \
  --target-msas data/target_msas.json \
  --output outputs/af3_pilot
```

An explicit `--single-sequence` alternative prepares empty target MSAs as well. Use that only as a separately documented protocol shared across arms. It cannot override a locked protocol. `--samples-per-design-seed N` retains candidate sample indices below N, and is an input-selection aid rather than a guarantee of equal generation cost. Freeze any selection before comparing results and retain the complete attempt ledger.

The prepared bundle contains input JSON, immutable reference copies and a hashed manifest. Output directories must be empty. No inference runs during preparation. The AF3 model may predict a different target arrangement; the design reference remains the experimental starting context and the proposed design pose.

## Execute in the AF3 environment

Install the pinned AF3 source with its supported runtime and obtain authorized model parameters. The local project environment does not install AF3 or fetch its weights. First inspect the launch plan:

```bash
python scripts/run_af3.py \
  --manifest outputs/af3_pilot/manifest.json \
  --af3-python /path/to/af3/environment/bin/python \
  --model-dir /path/to/af3/parameters
```

Add `--execute` to run. The launcher verifies the pinned source, hashes parameter files, records the resolved environment and launches only prepared inputs. It disables the data pipeline because MSA inputs are already supplied. Use `--resume --execute` only with an identical plan and engine. Do not mix external server predictions with this local-run record.

## Score and interpret

```bash
python scripts/score_af3.py \
  --manifest outputs/af3_pilot/manifest.json \
  --predictions outputs/af3_pilot/predictions \
  --output outputs/af3_pilot_scores
```

For each candidate, select the highest **unrounded AF3 ranking score**, then measure geometry on that selected sample. Never select the sample with the most favorable desired pose. All planned seed/sample outputs must be present; changed sequences, MSAs, templates, malformed confidences and extra samples invalidate the candidate result.

| Draft gate | Meaning |
|---|---|
| Target–binder pair ipTM ≥ 0.8 | Model confidence in the predicted interface |
| Binder fold Cα RMSD ≤ 2 Å | Consistency with the candidate's imported binder fold after binder alignment |
| Binder pose Cα RMSD ≤ 3 Å | Consistency with the candidate pose after target alignment |
| No severe clash flag | AF3's reported clash assessment |

The legacy `reference_coverage` field is an input completeness check. Current reference imports require a complete, finite backbone and the scorer rejects predictions with missing Cα atoms as invalid output. Coverage therefore equals 1.0 for every scored import, leaving four discriminating structural gates.

Also report target-alignment RMSD, contact retention, directional interface PAE and binder atom pLDDT. None is a binding-affinity measurement. The geometry reference is the candidate complex, so substantial generator-induced changes to the experimental receptor must be assessed independently before nomination. No direct 37NT deformation gate is implemented in the scorer.

Missing or invalid outputs remain in the planned denominator. `pass_fraction_of_planned` counts completed passes divided by all imported candidates, with incomplete outcomes explicitly labeled. `complete_cohort_pass_rate` is withheld until all arms are complete. `execution_verified` records matching local launch provenance; it does not establish biological validity. No independent-family inference, equivalence conclusion or automatic experimental nomination is made.

The summary describes the imported cohort, not all generator attempts. Compute-normalized yield requires external run logs. Confidence on MERTK combined with low confidence on AXL/TYRO3 does not prove selectivity; perform actual binding and functional assays. Multivalent constructs, the Annexin A5–PROS1 fusion and natural-ligand complexes need separate predictions and interpretation rather than forcing them through this two-chain schema.

Reference: [official AF3 output definitions](https://github.com/google-deepmind/alphafold3/blob/main/docs/output.md).
