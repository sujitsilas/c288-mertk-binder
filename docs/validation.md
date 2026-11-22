# Implementation and evidence status

## Current round: November 12 evening

The [new run plan](next_experiments.md) supersedes the historical two-attempt/two-sequence future-run recommendation below. The revised proposal keeps 20 MPNN proposals per sequence-design stage, changes epitopes, and uses a PD-L1 BindCraft control before releasing MERTK work.

All **48 local tests** pass, including independent-trajectory counting, censored-run handling, confidence bounds, live-budget parsing, and the real crop-neighborhood regression. Ruff, shell syntax and `git diff --check` pass. Frozen remote input validation job **1573710** completed with exit **0:0**, verifying 38 input/workflow hashes and the pinned engine source/runtime marker. These are software/input validations; new generation jobs are still queued as of November 12 at 10:53 p.m. Pacific.

The three GPU allocations are **1573755** (PD-L1 control), **1573761** (Ig1), and **1573763** (FNIIIb), eligible for A100/H200/32 GB V100. The control acceptance/geometry gate must pass before either MERTK allocation starts. No accepted control, MERTK hit or AF3 assessment is claimed. The manifest retains the original submissions and reconciled GPU-feature amendment. Source/proposal changes already present before this session were preserved.

Updated November 10, 2026.

## Verified locally

- `uv sync --frozen --extra dev` installs the lightweight project environment.
- All 21 tests pass (with the analysis extra installed). Unit/integration fixtures cover candidate import, residue/chain checks, cohort validation, immutable AF3 bundles, launch planning, fold-versus-pose geometry, sample selection, invalid/missing outputs and common scoring, independent generator backbones, unequal design yields and mismatched target/epitope rejection.
- A real download of experimental MERTK structure 37NT was parsed and exported: chain A contains 392 complete observed protein backbones, author residues 93–484.
- BindCraft job preparation succeeded on that real exported construct, using the pinned existing upstream checkout. This produced settings and a launch script without executing the design engine.

Reproduce local checks from the repository root:

```bash
uv sync --frozen --extra dev
uv run --frozen --extra dev ruff check src scripts tests
uv run --frozen --extra dev pytest -q
python scripts/fetch_reference.py --output data/references/37NT_check
python scripts/prepare_bindcraft.py \
  --target-pdb data/references/37NT_check/target.pdb \
  --output outputs/bindcraft_check --designs 1
```

See [the Hoffman2 pilot record](hoffman2_pilot.md) for installed remote environments, submitted jobs, model-weight provenance and queue status. GPU generation and assessment remain separate from local software validation.

Choose empty directories for repeated preparation. Fetch the pinned BindCraft source first if it is absent. The local structure preparation check does not approve an epitope or demonstrate a successful design.

## Not established

No MERTK design campaign, LigandMPNN refinement, actual AF3 GPU prediction, binding assay, agonism assay, shedding-protection experiment, efferocytosis improvement, resolution benefit or tissue-regeneration result has been completed by this workflow. Published success on other targets is not a MERTK result.

The RFD4-Proteina and LigandMPNN source revisions and upstream interfaces were checked, but their engines have not been installed or run locally. Their inference remains external to this package. AXL/TYRO3 negative-target scoring and glycan, membrane, ternary-ligand and multivalency analysis are planned biological review steps rather than automated features.

The repository keeps the experimental reference identity separate from each candidate's generated pose. The current scorer measures consistency with the latter; it does not automatically verify that a generator preserved the experimental target conformation. That review is required before ordering proteins.

The [proposal](repository_archive.md#proposal-drafts) and [experimental plan](experiment_plan.md) distinguish measurements needed for binding, agonism, receptor preservation, uptake, resolution and regeneration.

## Proposal and figure review, November 11, 2026

The revised submission renders to exactly two pages with the required five sections, 11-point body text and one-inch margins. Its compact schematic fits within the design-workflow allocation. The separate illustrated structure/hotspot guide renders to five pages; all pages were visually inspected for readability, overlap and clipping. It includes three complete, attributed CC BY figures from the 2026 structural study and our native-coordinate PyMOL views. Source image hashes and licensing are recorded in [figure sources](repository_archive.md#proposal-drafts).

The [assessment](repository_archive.md#proposal-drafts) distinguishes proposed anchor patches from validated hotspots, records glycan/species/selectivity concerns, and corrects confidence-score interpretations from the class notes. These are structural assessments, not binding or functional results. The 21 local tests and Ruff checks on the figure/PDF scripts pass.

## Focused course proposal, November 12, 2026

The [current proposal](repository_archive.md#proposal-drafts) narrows the class project to one monovalent human MERTK binding probe, with GAS6/PROS1 compatibility as its main biological test. A bounded BindCraft campaign and up to four nominations replace the broader tool-comparison, fusion, valency, shedding, mouse and tissue aims. The [scope decision](course_project_scope.md) defines proposed new nomination rules separately from the historical locked pilot. Those new rules have not yet been fully implemented or executed.

The updated submission renders to exactly two pages with the same five rubric sections, 11-point body text and one-inch margins. Section-overflow checks pass, and both pages were visually inspected. The original workflow figure and its generator were updated; the generator passes Ruff. This revision changes documents and figure generation, not the scientific inference pipelines or prior results.

## Bounded RFD4 secondary arm, November 12, 2026

At the user's request, the focused biological question now uses primary BindCraft generation and a small RFD4 secondary arm under the same 24-A100-hour course cap. The first tranche requests at most seven GPU allocation-hours. See the [new run record](hoffman2_course_pilot.md) for submitted IDs and dependencies; submission is not completed inference.

New code prepares immutable inputs, asserts RFD4 conditioning counts, gates MERTK generation on an upstream PPI control and screens raw outputs in native context without overwriting a deformed target. All 25 local tests pass, including four synthetic gate tests. The real 37NT crop retains 199 residues, unchanged sequence/backbone and author numbering 286-484. Remote checks verified input/workflow hashes, upstream revisions and the earlier BindCraft GPU-ready marker before submission. Full annotated-CIF validation is a submitted CPU stage, not a result inferred from file preparation. No AF3 prediction or new biological success is claimed.

## Ig1 extension and real condition audit, November 12, 2026

The first course CPU input job, 1530041, passed both upstream annotated-CIF validators but failed our additional assertion that parser residue IDs must remain author numbers 286–484. Its GPU dependencies were canceled before execution. It produced no new designs. The corrected checker explicitly maps ordered native identities and coordinates to author numbers, then checks every conditioned backbone and hotspot atom. It does not weaken the upstream validator, change receptor chemistry, or restore distorted generated coordinates.

Real reruns passed: Patch A input job **1537152**, exit 0:0, 199 residues / 796 backbone / 15 hotspot atoms; Ig1 job **1537154**, exit 0:0, 193 residues / 772 backbone / 13 hotspot atoms. Maximum native backbone deviations after rigid alignment are approximately 0.00102 and 0.00100 Å, consistent with PDB rounding. The saved maps confirm the materialized target uses residue IDs 1–199 or 1–193 rather than original author IDs. Each upstream input/control validation reports zero errors and four reviewed warnings: expandable-segment sentinel component, implicit indexing, unspecified formal-charge entries and absent explicit `_struct_conn`. Known nonzero charges and native coordinates are preserved; these warnings are not evidence of a validated biological interface.

All **39 local tests** pass, including residue-renumbering, altered hotspot identity, reordered sequence, distorted backbone, altered condition masks and conservative Slurm validation-state parsing. Ruff and shell syntax checks pass. The Ig1 dossier reproduces exposure, glycan clearance, human/mouse correspondence and an explicitly approximate GAS6–AXL overlay. The two new PyMOL PNGs were visually inspected, and the portable session is saved locally and in Downloads. See the [Ig1/cooperativity amendment](ig_binders_and_cooperativity.md). No measured binding, thermodynamic coupling, shedding protection or activation follows from these checks.

## Completed pilot review and future-run limits

Both new campaigns completed with exit 0:0 throughout, but yielded zero accepted BindCraft designs and zero RFD4 MERTK outputs passing native-context screening. The [results review](course_pilot_results.md) records actual counts, per-output geometry, omitted Ig2 R210 contacts and N170 glycan overlaps. Native outputs and final monitor reports are preserved; AF3 remains unsubmitted.

The original BindCraft limits counted relaxed trajectories and accepted variants rather than attempts and sampled sequences. Future preparation now sets `num_seqs=2`, and the project launcher limits all admitted design-loop attempts. Five additional tests, including a simulated upstream entry point, bring the suite to **44 passing tests**. Fresh preparation against the real Ig1 reference confirmed two sampled sequences and a two-attempt manifest limit. Ruff and shell syntax checks pass. These changes have not been tested in a new GPU allocation, deployed over the completed campaign workflow or used to resubmit a job.
