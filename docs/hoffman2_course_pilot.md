# Bounded course run: November 12, 2026

**Status update:** both fresh campaigns completed, with no accepted BindCraft design and no RFD4 candidate passing native-context screening. The [completed-run review](course_pilot_results.md) records outcomes, actual attempt/sequence counts, crop/glycan clashes and a correction to future BindCraft limits. The original `patchA_01` input audit failed before GPU generation. Submission tables below are historical records, not claims that jobs remain active.

Campaign `course_20261112_patchA_01` adds the requested RFD4 secondary arm to the focused human MERTK binding/ligand-compatibility project. Its purpose is feasibility and candidate generation, not an adequately powered method benchmark. The earlier `pilot_20261110_01` results and locked protocol remain unchanged.

## Location and resources

Remote root: `/u/scratch/s/sujit009/c288-mertk-binder`; campaign: `runs/course_20261112_patchA_01`. The user previously authorized `pi_pallard`; actual submissions confirmed account `pi_pallard`, partition `campus`, QoS `campus24`. Existing runtimes and cached weights are reused. No weight downloads or dependency installations were submitted.

| Stage | Submitted job | Limit | Dependency |
|---|---:|---|---|
| Input validation | 1530041 | 20 min, 8 CPUs, 16 GB, CPU only | None |
| BindCraft generation | 1530043 | 4 h, 8 CPUs, 64 GB, one A100 | Successful input validation |
| Upstream RFD4 PPI control | 1530048 | 1 h, 8 CPUs, 64 GB, one A100 | Successful input validation |
| RFD4 MERTK generation | 1530052 | 2 h, 8 CPUs, 64 GB, one A100 | Successful RFD4 control, including raw-geometry QC |
| Native output-quality report | 1530056 | 20 min, 8 CPUs, 16 GB, CPU only | Both MERTK generator jobs terminate, whether successful or failed |

These GPU allocations total at most **seven A100-hours**, within the course-wide shared 24-hour generation/control/AF3 cap. Wall-time limits are not predicted runtimes. No AF3 inference job is submitted without usable candidates and native-context review. Slurm's recorded states and logs, rather than this submission table, determine what has actually run.

## Same MERTK design question, different generators

Both MERTK arms use the native experimental human 37NT **A286-A484** FNIII-pair crop: 199 residues, an 80-residue binder, and provisional FNIIIa anchors **I365/W367/A369**. No receptor residues were redesigned or replaced. The complete A93-A484 reference and original glycan-bearing mmCIF are retained for context checks. Protein-only generation cannot establish mammalian glycan accessibility, ligand compatibility, agonism or shedding protection.

BindCraft is capped at **two trajectories**, at most two nominated MPNN sequence variants per trajectory, and two accepted designs. Upstream filters and optimization schedule are preserved; all saved trajectories and rejected/intermediate files remain available. Actual seeds are recovered from native names/logs. No additional standalone sequence-refinement arm is run.

RFD4 uses target backbone internal-coordinate conditioning, the full target sequence and all 15 finite sidechain heavy atoms across the three anchor residues as hotspots. Its MERTK stage requests **two replicates**, seed **202**, batch size one, `medium/front_loaded` recycling, no inline refolding oracle, and no compilation warm-up. Its Hydra snapshot is written inside the campaign. This is a distinct, cropped/hotspot-conditioned experiment, not a rerun of the earlier 392-residue unreviewed-surface input.

The RFD4 control is the pinned upstream PD-L1 PPI condition-spec example: target A17-A131, its original 50-120-residue binder range and hotspot atoms, one replicate, seed 101, batch size one, and the same medium/front-loaded schedule. Only the input-file path changes. The control tests basic runtime/conditioning geometry, not MERTK performance or the biological validity of an interface.

## Checks before and after generation

The prepared crop was checked against the original complete observed target: exact sequence, author numbering 286-484 and unchanged backbone coordinates. The upstream control contains 115 target residues, A17-A131. Original-reference hashes, staged input/workflow hashes and the pinned BindCraft/RFD4 source revisions were verified remotely before submission. BindCraft's prior GPU-ready record is required at launch.

The CPU input job materializes and checks both RFD4 annotated CIFs with the upstream inspector and full pipeline validator. A separate assertion checks 199 target residues, 796 conditioned N/CA/C/O atoms, complete target sequence conditioning and the intended 15 hotspot sidechain atoms. The general validator's reviewed chemistry warnings remain distinct from this condition-count assertion; neither test validates a therapeutic epitope.

The control must produce one complete two-chain structure with the exact target sequence, target C-alpha RMSD at most 2 Å after rigid alignment, no adjacent target/binder C-alpha separations outside 2.7-4.5 Å, and an inter-chain C-alpha contact within 8 Å. These are conservative diagnostic gates, not universal stereochemical criteria. A generated but geometry-rejected control exits with code 42 and prevents the RFD4 MERTK dependency from starting. This differs from a runtime crash; inspect `rfd4_control/quality.json` and the log.

The final CPU report screens native MERTK outputs against the same target/adjacent-C-alpha checks, intended-anchor contact and restored **context**: severe clashes with the complete protein or observed glycans are rejected. The raw generated target is never overwritten to manufacture a pass. Unresolved glycans and complete binder footprints still require review. Files/sequence variants are not counted as independent successful designs. Passing this report is not AF3 nomination, measured binding or biological success.

The RFD4 export code preserves a separate conditioned-coordinate reference annotation; that annotation is not proof that the actual generated coordinates obey the reference. The earlier raw files remain rejected. The new control is needed because the cause of their severe deformation remains unresolved.

## Monitoring and next step

Local manifest: `outputs/hpc/course_20261112_patchA_01/manifest.json`. It records source commit `8d1b141`, source revisions, input/workflow hashes, submission intents, exact Slurm commands and job IDs. Submission intents are written before requesting a job, so an uncertain SSH response must be reconciled rather than blindly resubmitted.

The Mac watcher polls every 120 seconds, saves scheduler states in `status.json`, and synchronizes logs, provenance and generated outputs. Once tracked jobs terminate, the existing plotter renders whatever real structures are available. `native_quality.json` is a separate quality report; diagnostic plots do not imply successful candidates. The local watcher requires the Mac to remain awake and connected; queued Slurm jobs do not depend on it.

```bash
.venv/bin/python scripts/hpc/monitor_pilots.py outputs/hpc/course_20261112_patchA_01/manifest.json --watch
```

After usable native outputs exist, review the complete glycan footprint, independence of candidates and native acceptance status. Then nominate up to four across both arms for the new common AF3 assessment described in the [scope decision](course_project_scope.md). The old cohort preparer requires both generators and a full-length target; it must not be reused unchanged for this cropped-target campaign. No empty or known-invalid AF3 job is queued.

The new preparation/quality scripts pass Ruff and shell syntax checks; all **25 local tests** pass. Four new synthetic tests cover chain assignment, receptor distortion, collapsed binder backbones and an observed-glycan clash. These checks validate workflow behavior, not the inference models or experimental biology.

## Continuation: corrected audit and an Ig1 surface

The original input job **1530041** failed after both upstream validators passed: our extra audit incorrectly assumed parser residue IDs retained the experimental author numbering. The three GPU dependencies were canceled without running; the final CPU quality report was empty. This was an input-audit error, not a failed binder-generation result. Original campaign files remain preserved.

Commit `90f2f45` replaces the assumption with an explicit map established by native ordered residue identities and backbone coordinates. It checks exact hotspot atom identities after mapping, plus unchanged sequence and coordinate conditioning. A separate submission check accepts Slurm's parsable output without requiring a trailing delimiter (`b657710`). No native chemistry or upstream validator was modified.

Two fresh campaigns use the same pinned runtimes, scratch root, account `pi_pallard`, partition `campus` and QoS `campus24`. Each CPU validation completed before its GPU jobs were submitted. Limits remain four hours BindCraft, one hour upstream RFD4 control and two hours RFD4 MERTK, one A100 / eight CPUs / 64 GB per GPU allocation. Each control gates its own MERTK RFD4 job. The combined GPU reservation is at most **14 A100-hours**, within the unchanged shared 24-hour course cap. The failed original tranche consumed no GPU allocation time.

| Stage | Patch A retry (`patchA_02`) | Ig1 side (`Ig1_01`) |
|---|---:|---:|
| Input validation | 1537152 | 1537154 |
| BindCraft generation | 1537698 | 1537774 |
| Upstream RFD4 PPI control | 1537710 | 1537775 |
| RFD4 MERTK generation | 1537712 | 1537780 |
| Native geometry report | 1537715 | 1537785 |

The Patch A audit verified 199 residues / 796 backbone atoms / 15 hotspot atoms; the Ig1 audit verified 193 / 772 / 13. Both exited 0:0. Parsed IDs are explicitly mapped back to native author numbers, and maximum aligned backbone deviations are approximately 0.001 Å. Both upstream validators reported zero errors per campaign input/control, with their four chemistry/indexing/sentinel warnings retained and reviewed. These are input checks, not predictions of binding.

The new Ig1 crop is A93-A285 with provisional K178/E183/V185 anchors. Patch A remains A286-A484 with I365/W367/A369. Both use 80-residue binders, two BindCraft trajectories and two RFD4 MERTK replicates. Full-receptor and observed-glycan context is screened using each campaign's own anchors. See the [structural dossier and co-binding/cooperativity plan](ig_binders_and_cooperativity.md). Sequence correspondence and an AXL ligand overlay are hypotheses; ligand compatibility, mouse binding and shedding protection remain unmeasured.

Independent Mac watchers poll and synchronize both campaigns every 120 seconds. They require an awake, connected Mac and do not repair or resubmit failures autonomously. Slurm jobs continue independently. The plotter runs after terminal job states and only plots available native structures. No AF3 inference is queued until real candidates pass review. As with the original table, these IDs record submissions; `status.json`, scheduler records and logs determine current execution state.

```bash
.venv/bin/python scripts/hpc/monitor_pilots.py outputs/hpc/course_20261112_patchA_02/manifest.json
.venv/bin/python scripts/hpc/monitor_pilots.py outputs/hpc/course_20261112_Ig1_01/manifest.json
```

Inputs and runtime workflow hashes were verified remotely before submission. Each manifest records base workflow commit `90f2f45`, submission-code identity, source revisions, exact scheduler commands and reconciled submission intents. All 39 local tests and Ruff checks pass. No new binder or biological success is inferred from these submissions.
