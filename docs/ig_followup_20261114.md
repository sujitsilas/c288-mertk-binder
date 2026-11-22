# Ig-only follow-up queued November 14, 2026

**November 16 evening update:** Both Ig generation jobs ran on V100 and reached their wall-time caps, preserving completed outputs. Ig1 produced two accepted variants from one backbone; both failed the frozen native-context screen despite clearing observed glycans. Ig2 produced no accepted variants. QC and protector triage completed, with no eligible parents. Parent AF3 job **1693273** was canceled while pending at **8:45 p.m. PST**, without allocating a GPU; the frozen workflow's empty-input branch recorded inference as skipped. The dependent CPU report then completed successfully at **8:47:56 p.m. PST**. No project job remains queued or running. TULP2 and mCherry AF3 had completed earlier on Blackwell, with all six planned outputs validated. The HTML now contains ten design poses and six comparator predictions; no parent is qualified.

Ig1 and Ig2 are distinct extracellular domains on every MERTK molecule. Binding both can produce **one-receptor occupancy**. An Ig1-Ig2 tandem is therefore a useful comparison, but cannot be called a receptor dimerizer merely because it has two arms. The alternative designs use two copies of an Ig1 binder or two copies of an Ig2 binder to test crosslinking of separate receptors. Linker geometry, preserved ligand function, receptor retention and actual apoptotic-cell internalization remain experimental questions. [TAM ectodomain structures](https://pmc.ncbi.nlm.nih.gov/articles/PMC13660871/).

The monovalent shedding-protector track shares the new Ig library. It does not duplicate generation with a different label. Its queued triage identifies provisional parents for subsequent assays of soluble MERTK, intact surface receptor and preserved efferocytosis. AF3 cannot establish cleavage inhibition. The published nanobody result supports testing monovalent protection separately from crosslinking; it does not establish protection by our chosen epitopes. [Duprez et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC13378718/).

**November 16 GPU-pool amendment:** The existing AF3 jobs 1693312 (mCherry), 1693284 (TULP2), and 1693273 (conditional parents) now also accept RTX PRO 6000 Blackwell nodes, retaining A100 80 GB and H200 eligibility. Their CPU request is four cores instead of eight to fit the available g36 CPU slots. Job IDs, queue age, memory, walltime, seeds and scientific inputs are preserved. No additional GPU allocation was submitted. The installed AF3 3.0.5 runtime uses JAX 0.10.2, Tokamax 0.0.12 and CUDA compiler 12.9; [official AF3 releases](https://github.com/google-deepmind/alphafold3/releases) document Blackwell support, and [CUDA 12.8 onward](https://docs.nvidia.com/cuda/archive/12.8.2/cuda-features-archive/index.html) compiles SM120. BindCraft's older compiler 12.6 remains outside that pool until separately validated. [Scheduler amendment](../outputs/maintenance/20261116/af3_gpu_amendment_20261116T204905Z/amendment.json).

**Account routing check:** `h2 whoami` confirms `pi_pallard` currently routes to campus because its group-owned nodes have not migrated to Slurm. Although the association lists `pi_shared24`, that alone does not establish access to the contributed-node pool. No shared-partition change was applied. Thus the currently spare L40S cards outside campus are not available to this account. [Official routing policy](https://www.hoffman2.idre.ucla.edu/Using-H2/Computing/Computing-Slurm.html).

## Submitted work

At 11:20 a.m. PST on November 14, input validation **1693261** had completed successfully in 10 seconds. The GPU preflight and comparator jobs were pending with `ReqNodeNotAvail`; downstream jobs were waiting on dependencies. Slurm gave no start-time estimate.

| Stage | Slurm job | Resource | Maximum allocation |
|---|---:|---|---|
| input_check | 1693261 | CPU | 00:10:00 |
| glycan_preflight | 1693262 | GPU | 00:20:00 |
| Ig1_E152_generate | 1693263 | GPU | 01:45:00 |
| Ig1_E152_quality | 1693265 | CPU | 00:15:00 |
| Ig2_L270_generate | 1693266 | GPU | 02:00:00 |
| Ig2_L270_quality | 1693268 | CPU | 00:15:00 |
| protector_triage | 1693269 | CPU | 00:10:00 |
| af3_parents | 1693273 | GPU | 00:45:00 |
| af3_tulp2 | 1693284 | GPU | 00:45:00 |
| af3_mcherry | 1693312 | GPU | 00:30:00 |
| af3_report | 1693319 | CPU | 00:15:00 |

The frozen suite is [round_20261114_ig_guard_01](../outputs/hpc/round_20261114_ig_guard_01/manifest.json). Submission accounting measured **17.5714 GPU-hours** spent/reserved before this suite. Its maximum was **6.0833 hours**, including two hours for AF3, for **23.6547/24 hours** total at submission. Following the preflight (104 seconds), TULP2 (195 seconds), mCherry (138 seconds), Ig1 (6,300 seconds) and Ig2 (7,202 seconds), tracked consumption is **21.4433/24 GPU-hours**. The empty parent-AF3 job used zero GPU time and has no remaining reservation; **2.5567 hours** remain under the project cap. These are allocated wall-clock GPU-hours; hardware throughput is not equated. BindCraft accepts A100/H200/V100; AF3 now accepts A100 80 GB, H200 141 GB or RTX PRO 6000 Blackwell.

Ig1 uses I122/Y148 and Ig2 uses L270, both on native A93-A285. Default BindCraft filters remain unchanged. The previous successful PD-L1 control is preserved as frozen evidence at the unchanged engine revision. The new generation objective requires a separate GPU integration preflight before either MERTK allocation may run.

## Measured generation results, November 16

| Run | Execution, PST | Completed attempts | Interrupted attempts | BindCraft-accepted variants | Native-context passes |
|---|---|---:|---:|---:|---:|
| Ig1 / 1693263 / V100 g9 | 4:44:13–6:29:13 p.m. | 5 | 1 | 2 from one backbone | 0 |
| Ig2 / 1693266 / V100 g8 | 5:04:11–7:04:13 p.m. | 9 | 1 | 0 | 0 |

Slurm reports `TIMEOUT` for both runs because they exhausted their declared allocations. The logs show ongoing optimization followed by time-limit termination; these were not out-of-memory or Python crashes. Their completed outputs were retained and QC ran afterward. Ig2 generated 40 saved sequence proposals; 38 failed base AF2 filtering, and neither of the two variants with full metrics passed all final filters. These small, time-limited samples do not establish that the Ig2 domain is intrinsically unsuitable.

The new Ig1 backbone is `Ig1_E152_l125_s438876`. Its `mpnn8` and `mpnn12` variants have mean AF2 iPTM **0.80** and **0.82**, respectively, and mean complex pLDDT **92**. Whole-Ig-pair alignment gives observed-glycan clearance of **5.12 Å** and **3.56 Å**, respectively, but protein minimum distances of **0.71 Å** and **1.50 Å**. The `mpnn8` target also exceeds the frozen 2 Å whole-pair RMSD threshold (2.07 Å). Both therefore remain unqualified.

A separate diagnostic using only Ig1 for alignment also finds protein overlap: **0.80 Å** for `mpnn8` and **1.57 Å** for `mpnn12`, while observed-glycan distances remain **3.94 Å** and **3.71 Å**. The nearest overlapping receptor atoms are side chains at **K109/K112**, plus **T145** for `mpnn12`; protein outside the designed Ig pair is more than 45 Å away. Thus the remaining rejection is not an omitted FNIII-domain collision. Whether native side-chain repacking can resolve these overlaps while retaining the interface and glycan clearance is untested. The original screening rules and qualification remain unchanged. This outcome does not isolate the causal contribution of the restraint from other changed settings. [Diagnostic atom-level audit](../outputs/maintenance/20261116/run_outcome_20261117T044539Z/Ig1_native_clash_audit.json).

Protector triage returned an empty parent list; no shedding-protective activity or receptor-dimerizing scaffold was established. The parent-AF3 skip is recorded in [execution_skipped.json](../outputs/hpc/round_20261114_ig_guard_01/assessment/execution_skipped.json). No new GPU experiment was submitted during this status check.

## Experimental glycan restraint

The adapter adds a differentiable excluded-volume loss to the pinned BindCraft hallucination stage. It transports **239 observed sugar heavy atoms at N114/N170/N207/N215/N234** using each Asn's N/CA/C frame, penalizes binder-heavy-atom proximity below **3.5 Å**, and applies a conservative **15 Å binder-CA/sequon-CA** avoidance term. The two weights are 10 and 2; the worst fifth of binder residues is averaged to reduce dilution with length. These are declared pilot settings, not calibrated physical energies.

Local tests cover collision versus distant controls, absent-atom masking, rigid-body invariance and native attachment mapping. The installed cluster JAX runtime passed finite-gradient, repulsive-gradient and zero-distant-penalty checks. **The one-step AF2 forward/backward GPU integration test passed on V100 g8 on November 14; job 1693262 completed at 9:29 p.m. PST in 104 seconds.** This verifies runtime integration, not biological glycan compatibility. A failed integration test blocks both generation jobs. The original upstream Python source is unchanged; the new wrapper and all constraint inputs are frozen and hashed.

This models a steric restraint around one observed insect-cell glycoform, not a complete mammalian glycan ensemble or joint glycan folding. Existing intact-receptor/glycan postchecks remain required. New binders may still fail them.

## Matched AF3 comparisons

TULP2 and mCherry are separate jobs with three seeds (101/202/303), one sample per seed, ten recycles, no templates and single-sequence MSAs. The target is the same observed human MERTK A93-A484 sequence. Glycans are not cofolded by this AF3 protocol. TULP2 is a reported Mer-Fc-negative comparator, while mCherry is an unrelated protein whose nonbinding in this system is unproven.

Both jobs completed on November 16: TULP2 used 195 seconds and mCherry used 138 seconds, each on one Blackwell GPU. Across their three seeds, TULP2 iPTM is **0.11–0.12** and partner atom-mean pLDDT is **29.45–30.89**; mCherry iPTM is **0.09–0.10** and partner atom-mean pLDDT is **31.55–33.33**. The interfaces and partner structures have low confidence under this single-sequence protocol. These results establish working inference and comparator outputs, not validated nonbinding or a calibrated binding assay. [Validated results manifest](../outputs/hpc/round_20261114_ig_guard_01/comparisons/results_manifest.json).

The parent assessment is conditional: at most one preliminary native-context-passing parent per Ig epitope, ranked by the recorded BindCraft interface pTM, enters the identical AF3 protocol. If no parent passes, inference is skipped and that absence is recorded. This is computational triage, not a native-binding nomination. All three outputs must be present with verified sequences and inputs before a case is displayed as complete. Low predicted confidence does not establish nonbinding.

The [Ig pilot gallery](prediction_review.md#rebuild-the-pilot-galleries) retains the eight rejected historical Ig1 poses and adds real new Ig outputs and comparator predictions as they synchronize. No placeholder protein pose is presented as a result. The local watcher runs every two minutes while the Mac is awake and connected; Slurm dependencies continue on the cluster independently.

## Repository cleanup

Unused local ProteinMPNN and RFD4 source checkouts and development caches were moved outside this repository to a recoverable archive. The remote runtimes remain installed. All frozen campaigns, measured failures, original structures, proposal artifacts and informative reports are retained. [Exact move list, hashes and recovery paths](../outputs/maintenance/20261114/cleanup.json).
