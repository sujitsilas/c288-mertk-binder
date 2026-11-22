# Workflow entry points

Run these commands from the repository root with the project environment installed.

| Task | Entry point |
| --- | --- |
| Fetch or inspect receptor references | [fetch_reference.py](fetch_reference.py), [review_ig_epitope.py](review_ig_epitope.py) |
| Prepare Aim 1 inputs | [prepare_ig_design_inputs.py](prepare_ig_design_inputs.py), [prepare_bindcraft_round.py](hpc/prepare_bindcraft_round.py) |
| Prepare Aim 2 inputs | [prepare_cleavage_protector.py](prepare_cleavage_protector.py) |
| Submit prepared Hoffman2 runs | [Ig follow-up](hpc/submit_ig_followup.py), [cleavage protector](hpc/submit_protector_round.py), [RFD4 repair](hpc/submit_rfd4_repair.py) |
| Run and score Boltz | [run_boltz_screen.sh](hpc/run_boltz_screen.sh), [score_boltz_screen.py](hpc/score_boltz_screen.py) |
| Analyze local Ig1 repairs | [repair_ig1_clash.py](hpc/repair_ig1_clash.py), [optimize_ig1_interface.py](hpc/optimize_ig1_interface.py) |
| Build the final HTML review | [build_final_review.py](build_final_review.py) |

The submission tools accept a frozen suite manifest and provide a dry run before submission. On Hoffman2, `squeue --me` shows live jobs and `sacct` reports elapsed time and exit status for the job IDs in a run. See the [compute guide](../docs/hoffman2_pilot.md) for historical setup commands.

## Regenerate viewers

```bash
.venv/bin/python scripts/build_final_review.py
```

For the standalone GitHub Pages export, run:

```bash
.venv/bin/python scripts/build_final_review.py --public --output output/site/index.html
```

The Pages source is the `gh-pages` branch. Its root contains the public HTML and `.nojekyll`. The build requires synchronized structure files. It embeds the Hoffman2 models and writes to `outputs/review/mertk_protector_predictions.html`. See [viewer provenance](../docs/prediction_review.md) for the original control inputs and on-demand pilot gallery builds.

## Analyze protease expression

Place the supplied H5AD at `data/burn_sham/filtered_final_burn_sham.h5ad`. The scripts also accept an explicit input path.

```bash
Rscript scripts/burn_sham/pseudobulk.R
Rscript scripts/burn_sham/pooled_transcripts.R
Rscript scripts/burn_sham/figures.R
```

The pseudobulk script uses sample-level summaries of normalized expression. The pooled analysis reconstructs transcript counts to measure abundance across the supplied cohort. Details are in the [protease analysis notes](../docs/burn_sham_mask_20261118.md).
