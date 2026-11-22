# Current queue and historical experiment record

**November 17 current work:** Independent cleavage-protector design jobs 1896138 and 1896184 are submitted. This creates an explicit stalk-targeting exception for Aim 2; the old F472 patch remains retired. [Triage and queue](cleavage_protector_20261117.md). Older dated updates below are historical.

**November 14:** The new Ig-only suite is submitted, including glycan preflight, Ig1, Ig2, protector triage and three AF3 jobs. [Current job IDs, dependencies, budget and method](ig_followup_20261114.md). The old Ig1 run finished with 23 complete attempts, one interrupted attempt and eight accepted variants from four backbones; all fail native glycan review.

The dated notes below are historical snapshots and do not describe the current queue.

# November 12 review and next experiments

**Prediction display and controls:** The [Ig pilot gallery](prediction_review.md#rebuild-the-pilot-galleries) displays synchronized Ig1/Ig2 poses with native receptor/glycan overlays and separate original-complex views. Following the user's request for different comparators, **TULP2** replaces Gokul's TULP3, and **mCherry** replaces albumin. Both control AF3 inputs are prepared, not queued. The current viewer has six filter-passing variants from three backbones; all six show native-context clashes. See the [control protocol and limits](prediction_review.md). The local watcher now refreshes the report after synchronization without submitting jobs.

**Latest clarification:** Two independent design aims: an Ig1/Ig2 MERTK binder and a MERTK cleavage protector. The current Ig1 baseline remains active; the latest synchronized metrics contain 14 starts, 13 completed attempts and four filter-passing variants from two backbones, all failing the N114 glycan review. [Two-aim plan](efferocytosis_design.md) · [Interim evidence](../outputs/analysis/20261113/ig1_two_aims_update.json).

**Active scope, November 13 evening:** FNIIIb is retired at the user's request. New MERTK epitopes are restricted to Ig1/Ig2, including the shedding-protection arm. The scheduler audit at 6:33–6:34 p.m. found no active FNIIIb jobs: 1573763 had timed out and 1573765 completed. Only Ig1 1573761 and its waiting QC 1573762 remained. No cancellation was necessary. Historical results and frozen inputs are retained; the old non-Ig submission routes are disabled in the current project launchers.

The next target bundle prepares **Ig1 I122/Y148** and **Ig2 L270**, preserving the full Ig pair, all 14 sequon locations and resolved glycan coordinates. L196 is masked because it is within 15 Å of N170. The [functional plan](efferocytosis_design.md) specifies M1/M2/M1–PS/M2–PS and a separate monovalent Ig shedding-protector screen. [Prepared inputs](../outputs/design/ig_only_20261113/manifest.json) are not GPU inference results. During-generation glycan exclusion remains an implementation prerequisite; no new Ig2, scaffold or protector jobs have been submitted. The four-hour proposed suite in the active config is subject to live accounting and constraint validation.

**Historical record below:** the original three-experiment round and earlier biological scope are retained for provenance; they do not authorize renewed FNIIIb sampling or direct stalk targeting.

**November 13 objective correction:** The [functional design plan](efferocytosis_design.md) makes improved MERTK-dependent efferocytosis the project endpoint. The frozen round below remains a binding-module discovery experiment; its filters and accounting are unchanged. An accepted binder, neutral probe or antagonist alone will not satisfy the functional objective.

**November 13, 6:11 p.m. update:** Ig1 **1573761** is running on a V100 32 GB (g9), with nine starts, eight completed attempts and one in progress in the synchronized snapshot. Two variants from one backbone pass BindCraft filters, but both overlap N114 carbohydrate and fail interim native-context review. No AF3 candidate is nominated. The [Abedi meeting follow-up](glycan_and_cleavage_followup.md) adds during-design glycan avoidance as a required next improvement and a separate monovalent cleavage-protector hypothesis. No new job or running-input change was made.

The revised proposal asks for selective MERTK binding and then classifies ligand competition; blocking is no longer a project failure. This record translates that change into three bounded BindCraft experiments on Hoffman2. The user subsequently requested GPUs beyond A100; the queued allocations were broadened to A100, H200 and 32 GB V100 without adding duplicate experiments or increasing the shared cap.

## FNIIIb diagnosis follow-up

The [structural diagnosis](fniiib_failure_diagnosis.md) identifies two independent problems in the only full-scoring FNIIIb family: all 13 variants fail the buried-polar cutoff, and their large interface overlaps the resolved N294 glycan even after fitting FNIIIb alone. The 37 parent interface positions were fixed across these variants. Prioritize the distinct Ig1 I122/F147/Y148 site already queued; this result does not establish that all FNIIIb epitopes are unsuitable. The earlier annotation of only four modeled glycan sites was also corrected to ten, consistent with the coordinates and structural paper.

## Status update: November 13, 8:04 a.m. Pacific

The PD-L1 control **1573755** completed in **10m42s** on an **A100 40 GB** and produced **two accepted variants from one independent trajectory**. Both passed the basic geometry check; control gate **1573758** completed successfully and released both MERTK dependencies. This establishes that the installed stack can produce accepted control designs; it does not establish a MERTK success rate.

FNIIIb **1573763** ran on the same GPU model until its four-hour allocation limit (**TIMEOUT**, elapsed 4h00m05s). It started **31 attempts**, completed **30**, and was interrupted during attempt 31. Fourteen trajectories reached structural scoring and produced **280 saved MPNN proposals**. **267** failed base AF2 validation; **13** reached full scoring but none passed the final filters. The overlapping failure counters include 267 interface-PAE failures, 219 interface-pTM failures and 13 interface-unsatisfied-H-bond failures. There are **zero accepted FNIIIb candidates**, so its completed native-context report contains no accepted structures to screen. An empty report is not a geometry success.

Ig1 **1573761** remains **PENDING** with `ReqNodeNotAvail`; its control dependency has cleared. It remains eligible for **A100, H200 and 32 GB V100**. Its quality job **1573762** waits for generation. No AF3 job has been submitted. The project has consumed **5.563 GPU allocation-hours** including historical runs, the new control and FNIIIb. Queued Ig1 time has not been consumed.

The latest logs, metrics and scheduler state were synchronized to `outputs/hpc/round_20261112_01/`. No jobs were resubmitted or resource limits changed during this status check. The local watcher remains alive; intermittent SSH/sleep interruptions have delayed some synchronizations, while Slurm execution continued independently.

## Evidence from completed experiments

| Metric | Historical Patch A | Historical Ig1 side |
|---|---:|---:|
| Trajectories actually started/completed | 4 / 4 | 6 / 6 |
| Trajectories with relaxed structural statistics | 2 | 2 |
| Saved MPNN sequence proposals | 40 | 40 |
| Proposals rejected at base AF2 validation | 40 | 40 |
| Accepted variants / independent trajectories | 0 / 0 | 0 / 0 |
| Interface-pTM failure counter | 40 | 34 |
| Interface-PAE failure counter | 40 | 38 |
| pLDDT failure counter | 15 | 36 |
| RFD4 native-context passes | 0 / 2 | 0 / 2 |
| One-sided exact 95% upper acceptance bound, per site | 52.7% | 39.3% |

Counts come from native generation logs, CSVs, saved FASTAs and Accepted directories, not from structure-gallery file counts. Failure categories overlap. The intervals concern independent design attempts under a fixed setting, not sampled sequence variants; these small samples do not constrain acceptance tightly. Do not pool different epitopes into one calibrated success rate.

Initial Ig1 trajectories had interface pTM 0.81/0.76, but no redesigned sequence passed base AF2 validation. Their surface hydrophobicities were 0.57 and 0.32; only the former exceeds the later 0.35 default threshold. This does not establish that missing aromatic anchors caused failure. Patch A's corresponding interface pTM values were 0.39/0.72 and hydrophobicities 0.40/0.38. Native trajectory pLDDT and PAE values use upstream normalized conventions; do not label them as AF3 pLDDT percentages or angstrom-valued interface PAE.

The six course GPU jobs consumed 4,427 seconds, or **1.230 A100-hours**. The accounting audit also includes 501 seconds of earlier RFD4 generation and 53 seconds of BindCraft GPU validation: **4,981 seconds = 1.384 hours** before this round. CPU setup time is excluded. Slurm elapsed allocation time is not GPU utilization or billing.

Run `scripts/summarize_campaigns.py` against the local campaigns to rebuild the metrics figure, JSON and CSV. The script requires an available generation log to count attempts. Supply `--log FILE` for a single campaign when its native log is missing and you have a specific scheduler log. The JSON retains source hashes. [Historical RFD4 results](course_pilot_results.md)

## Proposal review

The two-page PDF matches the revised Markdown and renders legibly. Both were retained as supplied; these are review notes rather than a silent rewrite of the submission.

- **Sampling is a goal bounded by resources.** At the previous end-to-end mean of 5.87 minutes per BindCraft start, 200 starts alone would need about 19.6 hours. Variable lengths, more successful sequences, different hardware and initial-guess validation can change that substantially. The primary allocation has a 12-hour cap, so 200 cannot be promised. Report actual complete and censored attempts.
- **A small positive control checks feasibility, not a published rate.** A matching accepted PD-L1 structure shows that this installed pipeline can reach its filters. Zero control hits in 20 attempts would still allow a one-sided 95% upper acceptance probability of 13.9%; it does not prove a broken installation. The gate pauses MERTK spending for diagnosis rather than relaxing filters.
- **The ranking is a heuristic.** Atomic-neighbor density is not a validated concavity metric. Glycan spheres approximate risk and are not explicit mammalian glycans. Re-running the script confirms 205 candidate patches, 32 spatial representatives and 14 observed-sequence sequons. The output now includes every patch and source hashes.
- **Neither new footprint is proven glycan-clear.** The Q124-centred patch has 31% of member C-alpha atoms within a 15 A sequon sphere; the V469-centred FNIIIb patch has 40%. Their patch centroids lie 15.8 A from N114 and 15.2 A from N442, respectively. The minimum member-heavy-atom distances to observed sugars are 8.71 A and 7.44 A. These are different geometric quantities. A patch centred on F472 itself is more shielded (89%), so the named FNIIIb arm must not be described as uniformly accessible.
- **Use the actual anchors.** The Ig1 job conditions I122/F147/Y148, the three exposed anchors selected by the screen; Y123 is not in this job's hotspot list. The FNIIIb job conditions V469/P471/F472. The 6.2 A and 64.4 A ligand-proximity figures refer to the respective patch centre residues and aligned functional positions, not a solved human ligand footprint. Antagonism and compatibility remain hypotheses.
- **Interpret a zero-hit bound precisely.** Zero successes in 200 completed fixed-setting trials gives a one-sided exact 95% upper bound of 1.49%. A two-sided interval is a different quantity. Early stopping, interrupted last attempts and adaptive conditions need separate accounting; the scripts report them explicitly.
- **Biological labels need defined uncertainty.** A competition ratio with a lower bound >=0.80 for both ligands supports the declared non-blocking criterion. For an antagonist label, require an upper bound <0.50 for the ligand in question and report each ligand separately. Intermediate or uncertain results remain unclassified. Failure on cells can implicate glycosylation, receptor presentation or assay differences; it does not identify glycosylation alone.

The [pinned BindCraft documentation](https://github.com/martinpacesa/BindCraft/blob/7713aa0d0d351e4117a8befeb8541f3a8ebd3368/README.md) supports initial-guess prediction for post-MPNN failures; it explicitly introduces pose bias. The [upstream wiki](https://github.com/martinpacesa/BindCraft/wiki/De-novo-binder-design-with-BindCraft) supports wider length sampling and warns of target-dependent sampling requirements. Neither source validates these MERTK sites. The [37NT record](https://www.rcsb.org/structure/37NT) identifies an insect-cell expression system and 3.41 A diffraction resolution, reinforcing the need for glycoform and local-interface caution.

## Submitted experiment specification

Suite: `round_20261112_01`. Scratch root: `/u/scratch/s/sujit009/c288-mertk-binder/runs/round_20261112_01`. Account `pi_pallard`, partition `campus`, QoS `campus24`; one GPU, eight CPUs and 64 GB host RAM per generation allocation.

| Experiment | Input / hotspots | Attempt cap | Binder lengths | GPU walltime | Job |
|---|---|---:|---|---|---:|
| PD-L1 positive control | Pinned upstream example; hotspot 56 | 20, or 2 accepted variants | 65-150 | 3 h | 1573755 |
| Ig1 primary | 37NT A93-A285; I122/F147/Y148 | 200 | 65-130 | 12 h | 1573761 |
| FNIIIb secondary | 37NT A286-A484; V469/P471/F472 | 40 | 65-130 | 4 h | 1573763 |

CPU input validation **1573710** completed before GPU submission. Control gate **1573758** follows control termination and exits nonzero unless an Accepted PDB has a matching native final-statistics row and passes the existing sequence/backbone/interface geometry check. Both MERTK generations depend on this gate succeeding. Their CPU native-context reports are **1573762** and **1573765**, respectively. Submission is not completed generation. At the last check (November 12, 10:55 p.m. Pacific), all generation jobs were pending: the control reported `ReqNodeNotAvail` and the MERTK jobs were waiting on dependencies. Slurm had no start-time estimate.

Every experiment samples 20 MPNN sequences per sequence-design stage and retains at most two accepted variants per trajectory. The PD-L1 control uses the upstream default multimer preset. Both MERTK arms use the corresponding hard-target preset with `predict_initial_guess=true`; this is a changed protocol, not an unbiased comparison to the old pilots. Default filters and optimization iteration counts are unchanged. Sampler seeds are 11012026, 11022026 and 11032026; actual sampled trajectory seeds remain in native names/logs. Exact deterministic replication across GPU models is not assumed.

The two MERTK crops preserve exact sequence, numbering and backbone coordinates. Their closest omitted receptor heavy atoms are **72.35 A** and **44.77 A** from their hotspot centroids, so both pass the 20 A context rule. This is a local input check, not proof that every eventual binder fits the complete receptor.

The maximum new reservation is **19 GPU-hours**. Including 1.384 hours already consumed and a **2-hour AF3 reserve** gives **22.384 of 24 hours**, leaving 1.616 hours uncommitted. All allowed GPU models count equally by elapsed allocation time for this cap; throughput and billing are not equated. No new RFD4 runs or AF3 inference were queued in this round. No project script automatically resubmits failed jobs.

The live campus inventory includes V100 32 GB, A100 40/80 GB and H200 141 GB cards. All are eligible through `A100|H200|V100`, following the user's broader-GPU request. The 8/11 GB cards are excluded for memory, and Blackwell RTX PRO 6000 cards are excluded from the pinned CUDA 12.6 stack pending separate compatibility validation. The [Hoffman2 Slurm guide](https://www.hoffman2.idre.ucla.edu/Using-H2/Computing/Computing-Slurm.html#requesting-gpus-with-slurm) documents model features and GPU requests. Available inventory does not guarantee immediate scheduling.

## Assessment and next decision

1. Read the control gate and its native metrics. If it produces zero accepted structures, inspect the failed stage and actual sample size before preparing a new, separately recorded experiment. MERTK jobs will be blocked/canceled by dependency rules.
2. For each MERTK run, count starts, complete attempts, censored last attempts, relaxed trajectories, saved sequence proposals, accepted variants and independent successful trajectories. Preserve overlapping failure counters and native output paths. The wrapper now consumes the manifest's attempt limit instead of the old hard-coded two-attempt limit.
3. Screen **accepted** MERTK complexes against complete A93-A484 context and resolved glycans. Review unresolved sequons, intended anchor recovery, domain preservation, crop-edge artifacts, exposed hydrophobicity and buried polar groups. The automatic screen is preliminary; a passing geometry score is not nomination.
4. Only then nominate up to six diverse, independent trajectories for an AF3 assessment on the complete observed receptor with seeds 101/202/303. Retain the separate declared binder-confidence, pair-ipTM, interface-PAE and pose-recovery gates. The old full-length/two-generator AF3 preparer does not implement this cropped, BindCraft-only nomination plan and must not be used unchanged. No AF3 job is automatically released by the new suite.
5. Characterized binding candidates, if any, enter the proposal's expression/SEC/stability, mammalian-MERTK binding, AXL/TYRO3 selectivity, cell-binding and ligand-competition sequence. Up to six is a ceiling, not a quota; no proteins have been nominated or ordered.

## Monitoring and reproducibility

```bash
.venv/bin/python scripts/hpc/monitor_pilots.py outputs/hpc/round_20261112_01/manifest.json --watch
```

The frozen suite contains 38 hashed input/workflow files, per-experiment manifests and a copy of the proposal at preparation. The original A100 submission commands are retained, with a separate reconciled scheduler amendment for the broadened GPU features. `gpu.txt` in each experiment records actual allocated hardware. Repo edits do not alter the frozen queued workflow. The local watcher was started and its process record is `monitor_process.json` beside the suite manifest. Local monitoring requires an awake, connected Mac; scheduler dependencies do not. The account permits at most two concurrently allocated GPUs, so a broader pool improves placement options without bypassing that limit.

Completed pilot data, native outputs, reference structures, weights and installed environments were retained. Cleanup moved 127 files (about 6.76 MiB) of duplicate log copies, obsolete render inspections, old preparation smoke checks and regenerable caches into the Mac Trash. The exact recovery paths and hashes are recorded in `outputs/maintenance/20261112/cleanup.json`; the user's pre-existing working-tree changes were also backed up there before edits. Monitoring now retrieves only the tracked jobs' logs to avoid recreating those duplicates.
