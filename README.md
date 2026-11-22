# MerTK binder design

This C288 project tests two approaches to preserving efferocytosis during inflammation. Aim 1 designs binders for human MerTK Ig1. Aim 2 tests whether a binder near the receptor stalk could protect it from shedding.

[Interactive results and 3D models](https://sujitsilas.github.io/c288-mertk-binder/) · [Design and experimental plan](docs/efferocytosis_design.md) · [Methods and run history](docs/README.md)

## Repository layout

- `src/mertk/` contains structural analysis and target-scope checks. `tests/` covers their use in the design workflows.
- `scripts/` contains structure preparation and rendering. `scripts/hpc/` contains Hoffman2 workflows. `scripts/burn_sham/` contains transcript analysis.
- `configs/` specifies targets and workflow settings. `assets/` contains the page template, figures, and pinned 3Dmol.js library.
- `data/`, `outputs/`, and `Report/` are local. They hold research inputs, synchronized run outputs, and the submission document. Run paths cited in historical notes refer to the local workspace.

## Reproduce the review

```bash
uv sync --frozen --extra dev
.venv/bin/pytest -q
.venv/bin/ruff check src scripts tests
.venv/bin/python scripts/build_final_review.py
```

Tests that use experimental structures run when those references are installed. The review build reads the synchronized runs and creates `outputs/review/mertk_protector_predictions.html`. From a clean checkout, restore and verify the bundled evidence before building. [Clean checkout instructions](docs/prediction_review.md#rebuild-the-published-review) · [Pilot gallery rebuilds](docs/prediction_review.md#rebuild-the-pilot-galleries) · [Archived material](docs/repository_archive.md)
