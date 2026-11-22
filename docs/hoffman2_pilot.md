# Hoffman2 pilot: 2026-11-10

The new [November 12 course campaign](hoffman2_course_pilot.md) uses a separate cropped-target Patch A pilot, a primary BindCraft arm, and a geometry-gated RFD4 secondary arm. This document retains the original campaign and its failures.

Campaign: `pilot_20261110_01`. The user authorized `pi_pallard` after checking the available Slurm accounts. `sujit009` is the login, not a registered Slurm allocation account. This campaign uses the production Slurm host `hoffman2.oarc.ucla.edu` through the existing SSH configuration; no passwords are stored in this repository.

## Storage and environments

Everything remote is under `/u/scratch/s/sujit009/c288-mertk-binder`:

| Content | Relative location |
|---|---|
| Pinned upstream source | `external/BindCraft`, `external/RFD4-Proteina`, `external/alphafold3` |
| BindCraft runtime | `envs/bindcraft-wheels` (recovery installation); Python bootstrap in `envs/miniforge/envs/BindCraft` |
| RFD4-Proteina environment | `external/RFD4-Proteina/.venv` |
| AF3 environment | `envs/alphafold3` |
| BindCraft AF2 weights | `external/BindCraft/params` |
| RFD4-Proteina release bundle | `cache/huggingface/hub/models--nvidia--RFD4-Proteina` |
| Standard AF3 parameters | `weights/alphafold3`, direct official Google download |
| Inputs and generated results | `runs/pilot_20261110_01` |
| Slurm logs | `logs` |
| Resolved package inventories | `provenance` |

The RFD4 release bundle is pinned to Hugging Face revision `a06e2509c7033b010fc43614e9fd840303588502`. Upstream code revisions are in `configs/external_tools.json`. The local ignored campaign manifest preserves inputs, commands, job IDs and retry history. Runtime caches and temporary directories are also assigned to scratch. No model weights are downloaded to the Mac.

## Pilot design and resource limits

Both generators use experimental 37NT chain A, author residues 93–484, and an 80-residue binder. This is an unreviewed surface-exploration/runtime pilot, not a selected therapeutic epitope. The protein-only input omits glycans and the unresolved membrane-proximal stalk; therefore this run cannot establish ligand compatibility, shedding protection, receptor clustering or agonism.

- BindCraft: at most two trajectories, stop after one accepted design, upstream default filters and four-stage schedule. Actual random seeds come from trajectory names/logs. Preserve rejected structures and logs.
- RFD4-Proteina: two replicates, seed 101, batch size one, `medium/front_loaded` recycling, target backbone internal-coordinate conditioning plus complete target sequence. No hotspot constraints, refolding oracle or shared confidence metrics. Input conversion and validation precede inference.
- AF3: the pinned runtime installed successfully. Current official instructions provide the standard parameters directly from Google; download job `1314141` completed successfully, storing them in scratch, validating the complete Zstandard stream, and recording the immutable object generation and SHA256 (`74d0258616917cd122f5eab6d076afe4a8930e96823851e65e4f777dfb1f33ff`). The earlier assumption that the user needed to supply access was corrected after checking the updated README. Assessment was initially queued after candidate selection from both generators, then canceled when preparation failed (see results below). The locked `configs/af3_pilot_single_sequence.json` uses identical single-sequence inputs, seed 101, five samples and ten recycles for both arms, with no structural templates. This is a runtime pilot; its lack of target MSA and tiny cohort preclude a broad method-performance claim.

Each setup requests 8 CPUs, 32 GB RAM and four hours with no GPU. Each generation job requests one A100, 8 CPUs, 64 GB RAM and two hours. Hoffman2 routes this account to `campus` with `campus24` QoS. Generation depends on successful setup and is canceled if that dependency fails. No full genetic database download is scheduled. A separate 15-minute, 4-CPU, 16-GB input-validation job checks the RFD4 condition file while the GPU pilot waits in the queue.

## November 11 result review

The first complete overnight attempt produced RFD4 files but no shared AF3
assessment. BindCraft setup `1315796` spent its four-hour allocation resolving
Conda requirements and timed out; generation `1315800` was canceled without
running. Candidate preparation `1315677` failed because BindCraft had no usable
complex, and dependent AF3 job `1315742` was canceled. AF3's runtime and weights
remain installed in scratch.

RFD4 job `1310910` completed in 8 minutes 21 seconds on an A100 and wrote both
requested 80-residue binder complexes. Completion is not design success:

| Native output | MERTK C-alpha RMSD to 37NT | Target residues with a binder C-alpha within 8 Å |
|---|---:|---:|
| Replicate 0 | 23.264 Å | 14 |
| Replicate 1 | 21.294 Å | 0 |

The input target matches 37NT to 0.00009 Å after rigid alignment (serialization
rounding). Both output files retain the full target sequence and annotations for
1,568 conditioned backbone atoms. The severe output deformation is real; its
mechanistic cause is unresolved. Neither file is a validated binder. The log
also records multiprocessing temporary-directory cleanup errors, but generation
finished, both output structures parse, and Slurm recorded exit code 0.

The original assessment policy selected replicate 0 by sequence/contact checks.
That selection and all original files are preserved. **After reviewing this
failure**, future candidate preparation additionally requires target C-alpha RMSD
at most 2 Å from experimental 37NT after rigid alignment. This is a newly added
structural quality gate, not a retrospectively preregistered success criterion;
both existing RFD4 outputs fail it. AF3 will not be requeued on this cohort.

BindCraft's recovery narrows the upstream dependency ranges to JAX/JAXlib 0.4.35,
CUDA 12.6, cuDNN 9.3, NumPy 1.26.4, Haiku 0.0.13, Flax 0.9.0, Optax 0.2.4 and
Chex 0.1.87. The specific conda-forge Python 3.10 GPU build of JAXlib requires
CUDA >=12.6; the observed A100 driver supports CUDA 13.2. The package-manager
step now has a 30-minute timeout inside the unchanged four-hour setup allocation.
Replacement setup job `1429715` and dependent generation job `1429716` were
submitted under `pi_pallard`, routed to `campus` / `campus24`, with the same
CPU, memory, GPU and wall-time limits. Local two-minute monitoring was restarted.
This is a bounded installation retry, not a verified repaired runtime. The
upstream source and scientific settings are unchanged. See the
[package metadata](https://api.anaconda.org/release/conda-forge/jaxlib/0.4.35)
and [upstream installer](https://github.com/martinpacesa/BindCraft/blob/7713aa0d0d351e4117a8befeb8541f3a8ebd3368/install_bindcraft.sh).

The Mac monitor recovered from intermittent SSH failures, synchronized the
terminal states and outputs, rendered three PNG/SVG diagnostic figures, then
exited as designed. Its status now describes AF3's observed scheduler state,
instead of repeating the original queued intent after cancellation. Figures are
under `outputs/hpc/pilot_20261110_01/figures/`; they show native predictions only.

## November 11 installation recovery and structure review

The narrowed Conda retry `1429715` reached its 30-minute package-manager limit;
its generation dependency `1429716` was canceled. The alternative CPU dependency
check `1459280` completed in 39 seconds. Recovery setup `1459533` runs under
`pi_pallard` with 8 CPUs, 32 GB and four hours, no GPU, in the same scratch
project. It uses an isolated Python 3.10 venv, upstream-compatible JAX 0.4.35 CUDA
wheels and Flax 0.9.0, with PyRosetta and NumPy 1.26.4 resolved together.
The initial wheel recovery `1459533` installed packages and passed dependency
checking but failed CPU import because the newer NVCC namespace package has
`__file__=None`, which JAX 0.4.35 passes to `Path`. The retry pins
`nvidia-cuda-nvcc-cu12==12.6.77` and reuses the environment.
ColabDesign is pinned to `8fd167defd1b0f98e29f00584ff4f3987e6935dc` and PyRosetta
to `2026.29+releasequarterly.80a0635615`. Source is unchanged. The installer
checks dependencies/imports, fetches the upstream AF2 parameters into scratch,
and records the environment. GPU validation remains a separate step; the ready
record explicitly sets `gpu_validated=false`. No replacement GPU generation was
submitted during surface review.

The [local PyMOL review](structure_review.md) includes experimental MERTK and
GAS6, predicted full PROS1/GAS6, experimental ligand-context fragments, two
exploratory surface patches and failed RFD4 overlays. It documents severe
within-domain distortions and adjacent C-alpha anomalies as well as whole-target
RMSD. The underlying RFD4 cause remains unresolved; do not infer a checkpoint or
conditioning bug solely from these measurements. Next diagnosis should compare
processed conditioning and reconstruction against an upstream PPI control
before using further MERTK generation outputs. User-facing coordinates and the
session are local; all model weights remain in scratch.

## November 12 status check

Remote Slurm accounting and the actual job log confirm that BindCraft GPU check
`1462338` **completed successfully**, exit `0:0`, on A100 node `g15` from
08:57:17 to 08:58:10 Pacific on November 12 (53 seconds). The ready record now
has `gpu_validated=true`. BindCraft function imports, JAX GPU matrix
multiplication/convolution and parameter-file presence checks passed. This is
a runtime smoke test; it did not run binder optimization or prove the complete
design workflow.

Recovery installation `1459732` completed in 3 minutes 49 seconds on November
11. The original Conda setup `1429715` remains failed, and its generation job
`1429716` was canceled before starting. No replacement generation job was
submitted during the structural review. There are still no BindCraft-generated
candidates in this campaign. RFD4's two native outputs remain rejected for
receptor deformation; AF3 assessment remains canceled, with no predictions.

The local monitor recorded all tracked jobs as terminal at 08:59 Pacific,
synchronized the successful check and regenerated the existing three diagnostic
figures before exiting normally. These figures still describe the two original
RFD4 structures, not new designs. A direct `squeue` query during this check timed
out at the scheduler, so no claim is made about unrelated jobs in the user's
live queue; `sacct` and remote log/provenance reads succeeded.

Next work is a bounded BindCraft generation run using a documented surface
specification and a separate RFD4 conditioning/reconstruction diagnostic. The
structural review has not validated either proposed patch as a production
hotspot. No jobs were resubmitted as part of this status check.

## Original jobs and observed setup issues

| Stage | Current job ID |
|---|---|
| BindCraft setup | `1315796` |
| BindCraft generation | `1315800` |
| RFD4-Proteina setup | `1310909` |
| RFD4-Proteina generation | `1310910` |
| AF3 environment setup | `1310912` (completed) |
| AF3 parameter download | `1314141` (completed) |
| Candidate selection / AF3 preparation | `1315677` |
| Shared AF3 assessment and scoring | `1315742` |
| RFD4 input validation (CPU) | `1313033` |

The first BindCraft setup (`1310891`) failed because the unversioned Miniforge installer alias lacked a checksum asset. Its GPU dependency (`1310908`) was canceled without running. The versioned installer and checksum verified successfully on retry. The retry generation submission response timed out, but Slurm accepted job `1311330`; its account, dependency and resource limits were verified with `scontrol` instead of submitting a duplicate.

BindCraft retry `1311229` created its initial Python environment under `envs/BindCraft`, still in scratch, because Mamba inferred a different root from the nested installation path. Its dependent generation `1311330` was canceled. The next attempt explicitly sets `MAMBA_ROOT_PREFIX` and `CONDA_ENVS_PATH` to the private Miniforge base. That attempt (`1313456`) was then stopped after prolonged dependency solving across inherited extra channels. The current setup (`1315796`) reuses the Python 3.10 environment and restricts Mamba to the upstream conda-forge/nvidia channels with strict channel priority; it leaves the user's global configuration unchanged.

The initial RFD4 input check (`1312273`) rejected PDB input because the condition-spec interface requires mmCIF chemistry metadata. The input now uses the original experimental `37NT.cif`, cropped to the same A93–484 protein segment; the GPU job depends on the corrected input check (`1313033`), which completed with zero errors. It confirmed 392 target sequence residues, 1,568 conditioned backbone atoms and an 80-residue expandable binder. Four warnings were reviewed: the intentional binder placeholder, implicit all-indexed default, zero formal charges serialized as `?` by the installed Biotite writer and restored as zero by its reader, and absence of polymer/non-polymer covalent link metadata in this protein-only input. The ignored `inputs/review.json` records this review and the remaining glycan/disulfide limitations. No charges were changed.

RFD4-Proteina environment installation, dependency checking and the pinned release-weight download completed successfully. GPU inference and remaining installations must be checked in the live status rather than inferred from this document.

## Monitoring and figures

From the Mac repository root:

```bash
uv sync --frozen --extra dev --extra analysis
python scripts/hpc/monitor_pilots.py outputs/hpc/pilot_20261110_01/manifest.json --watch
```

The monitor polls every two minutes, writes `status.json`, and copies only this project's logs, run outputs and environment inventories. It reloads the manifest each time so recorded replacement job IDs are followed. It does not retry failed scientific jobs automatically. It needs the Mac to remain awake with SSH/network access; the Slurm jobs continue independently.

When all recorded jobs are terminal, it runs `scripts/plot_pilot.py`. A fourth confidence-versus-pose figure is added only when verified shared AF3 scoring records exist. Real complexes produce PNG and SVG figures for a structure gallery, receptor contact profile and target-deformation/contact overview. No figures are fabricated if there are no generated complexes. Exact duplicate files are excluded, while skipped files are recorded. Intermediate/rejected designs and relaxed/sequence variants are explicitly distinguished from independent successful designs. All geometric measurements, hashes and chain assignments are saved in `figures/geometry.json`.

Twenty-one local tests pass, including synthetic checks of chain assignment, geometric calculations, duplicate/mismatched-target accounting the no-output case, and selection of one complex per BindCraft trajectory seed for a matched two-arm AF3 cohort. These tests validate analysis code, not model inference or biological activity.

## Preparing a fresh repeat

Use `scripts/hpc/prepare_pilot.py --campaign NEW_ID --remote-project /u/scratch/s/sujit009/c288-mertk-binder` after fetching the pinned source and experimental reference. It verifies reference hashes, sequence and numbering, copies the original mmCIF for RFD4 and protein PDB for BindCraft, and records the exact settings and input hashes. It refuses to overwrite an existing campaign. Stage those files to the recorded remote directory before submission. `scripts/hpc/submit_pilots.py MANIFEST --account ACCOUNT` submits actual jobs and requires an explicit account.

AF3 parameter source: [official acquisition instructions](https://github.com/google-deepmind/alphafold3#obtaining-model-parameters). The standard Google Storage object generation is `1780568696389861` (1,020,545,840 compressed bytes). Model parameters remain in scratch and are never synchronized to the Mac or committed.

## Dependent AF3 stage

The candidate-preparation job runs after both generators terminate, preserving any complete partial outputs. It selects at most two candidates per arm. BindCraft accepted complexes take precedence over relaxed/raw trajectory complexes, with at most one complex per actual native seed. RFD4 generated complexes are considered in sorted file order. Both arms must preserve the complete target sequence, contain complete backbones and have an inter-chain C-alpha contact within 8 Å. Selected and rejected paths, native acceptance status, hashes and seeds are recorded. This selected runtime cohort must not be used to estimate unbiased generator success rates.

AF3 is submitted with an `afterok` dependency on successful cohort preparation; a missing arm cancels this assessment rather than manufacturing a comparison. The same job scores prediction completeness and structural consistency, retaining missing/invalid predictions in the planned-candidate denominator. The Mac monitor downloads these reports and plots available data after all stages finish.

## Checking jobs yourself

```bash
ssh sujit009@hoffman2.oarc.ucla.edu
squeue --me
h2 why-pending 1310910
h2 diagnose 1310910
sacct -j 1310910 --format=JobID,JobName,State,Elapsed,ExitCode
```

`Priority` means the job is waiting its turn; `Dependency` means an earlier pipeline stage must finish. The original local `hoffman2` SSH alias still targets the legacy hostname, so use the new hostname explicitly. Source: [UCLA's Slurm guide](https://www.hoffman2.idre.ucla.edu/Using-H2/Computing/Computing-Slurm.html).
