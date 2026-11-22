# Patch A and Ig1 pilot results

> Historical run review. The [next round](next_experiments.md) restores 20-sequence sampling, uses new sites and a BindCraft positive control, and supersedes the two-attempt/two-sequence future-run recommendation below. Native pilot records remain unchanged.

Both November 12, 2026 campaigns completed. Every input, generation, control and quality-report Slurm allocation exited **0:0**. This establishes completed execution; **neither campaign produced a candidate passing the present selection criteria**. No AF3 prediction was submitted. Historical inputs and native outputs remain unchanged.

## What actually ran

| Measurement | Patch A, `course_20261112_patchA_02` | Ig1, `course_20261112_Ig1_01` |
|---|---:|---:|
| BindCraft job | 1537698, completed in 25m32s | 1537774, completed in 33m09s |
| Actual logged BindCraft trajectory starts | 4 | 6 |
| Trajectories reaching the recorded relaxation/sequence-design stage | 2 | 2 |
| Saved MPNN sequence proposals | 40 | 40 |
| Accepted BindCraft designs | **0** | **0** |
| RFD4 MERTK job | 1537712, completed in 1m45s | 1537780, completed in 9m34s |
| RFD4 MERTK replicates generated | 2 | 2 |
| RFD4 replicates passing the native-context screen | **0** | **0** |
| Upstream RFD4 PD-L1 control | Passed, 1537710, 2m08s | Passed, 1537775, 1m39s |

The six one-A100 allocations consumed a combined **4,427 GPU allocation-seconds = 1.230 A100-hours**, excluding CPU-only input/QC jobs and queue waiting. Their aggregate wall-time reservation had been 14 hours. These are elapsed allocation hours, not a utilization or billing measurement. The original `patchA_01` GPU jobs were canceled before running.

Counts distinguish attempts, completed trajectories, sequence proposals and saved structure files. BindCraft's final log says five/seven trajectories because its counter starts at one and increments after each loop; actual `Starting trajectory:` entries are four/six. The gallery contains eight Patch A and ten Ig1 complex **files**, including rejected and relaxed duplicates of the same trajectory; those are not eight/ten independent candidates.

## Why there are no accepted candidates

BindCraft reached two relaxed trajectories per epitope, but **all 40 MPNN sequence proposals per epitope failed its base AF2 validation filters**. Patch A's failure counters report 40 interface-pTM failures, 40 interface-PAE failures and 15 pLDDT failures. Ig1 reports 34, 38 and 36, respectively. These counts overlap and must not be added as independent failures. Confidence during initial backbone optimization is not the subsequent sequence-refolding result.

Both upstream RFD4 controls passed, with target C-alpha RMSD approximately **0.159 Å** and normal adjacent C-alpha distances. All four new MERTK replicates likewise pass the basic backbone geometry checks, including target C-alpha RMSD below 2 Å. The severe target deformation in the old full-length pilot is not reproduced in these cropped runs; this does not identify the earlier deformation's cause.

| RFD4 output | Target C-alpha RMSD, Å | Closest native-protein contact after alignment, Å | Closest observed-glycan contact, Å | Screen outcome |
|---|---:|---:|---:|---|
| Patch A rep0 | 1.386 | 1.328 | 4.364 | Protein overlap |
| Patch A rep1 | 1.002 | 0.213 | 1.541 | Protein and glycan overlaps |
| Ig1 rep0 | 0.757 | 2.110 | 0.490 | Glycan overlap |
| Ig1 rep1 | 0.897 | 0.553 | 0.568 | Protein and glycan overlaps |

Ig1 RMSDs above are rounded from the saved report; consult `native_quality.json` for exact values. Screening aligns the raw generated target to its experimental crop, then compares the transformed binder with the full experimental receptor and its observed glycans. It does not replace generated coordinates. Its conservative distance gates are protein minimum 2 Å and observed-glycan minimum 2.5 Å, together with backbone and intended-anchor criteria.

The two Patch A binders clash specifically with **Ig2 R210**, outside their FNIII-only design crop. Rep0's closest pair is binder T48 OG1 to native R210 NH1 (1.328 Å); rep1's is binder T58 O to R210 NE (0.213 Å). Thus the omitted Ig region matters even though the intended epitope lies in FNIIIa.

Both Ig1 binders overlap carbohydrate chain C, attached to **N170**: rep0 binder residue 1 O to NAG2 C2 (0.490 Å); rep1 binder residue 24 CG to NAG2 O6 (0.568 Å). The first NAG C1 lies 1.434 Å from N170 ND2, supporting that attachment assignment. Anchor-only glycan clearance did not guarantee space for an entire binder. The Ig1 BindCraft trajectory files also fail the intended-anchor/contact-context review; initial hotspot conditioning is a soft design objective rather than a demonstrated recovered interface.

These are **native-context screening failures**, not uniformly internally collapsed complexes. In their own generated receptor coordinates, the four RFD4 structures have minimum interchain heavy-atom distances of 3.280, 3.030, 2.263 and 2.293 Å, respectively. Restoring experimental sidechains/domain context can introduce contacts absent in the predicted complex. Preserve that distinction and review geometry locally before deciding whether a future design is salvageable. No passing candidate is declared by changing alignment or screening thresholds after seeing these results.

## A budgeting mistake and its correction

Our initial preparation used upstream `max_trajectories=2` and `max_mpnn_sequences=2`. Inspection of the pinned implementation shows these limit **relaxed trajectories** and **accepted** MPNN variants, respectively. They do not cap attempted trajectories or sampled sequences. `num_seqs` remained 20. Consequently, the promised two-attempt/two-sampled-sequence limits were not met: actual attempts were four/six and sampled sequences forty/forty. The scheduler time limits held, and these counts remain in the results.

For **future campaigns**, preparation now explicitly sets `num_seqs=2`. The project launcher `run_bounded_bindcraft.py` wraps the upstream loop's stopping check in memory to count every admitted loop iteration, including rejected trials, with a two-attempt maximum. It also respects upstream stopping conditions and refuses a nonempty design directory. The pinned upstream checkout, optimization schedule and filters are unchanged. Existing campaigns are not edited or silently rerun. A simulated upstream entry-point test checks that the real launcher path enforces two failed attempts; this is not a new GPU validation.

## Next decision and available figures

Review a full-receptor or appropriately expanded context for Patch A and the complete N170 glycan footprint for Ig1 before preparing another bounded batch. Simply increasing the number of the same cropped/protein-only inputs risks repeating these failures. Anchor exposure and a small unsuccessful pilot do not establish that either epitope is intrinsically inaccessible or that one generator is superior.

The Mac monitors synchronized final files and generated three PNG/SVG figures per campaign under:

- `outputs/hpc/course_20261112_patchA_02/figures/`
- `outputs/hpc/course_20261112_Ig1_01/figures/`

`complex_gallery`, `target_contact_profile` and `geometry_overview` are inventories and geometry diagnostics, not rankings or confidence/affinity figures. Their `geometry.json` files retain hashes and skipped files. Detailed selection reasons are in each campaign's `native_quality.json`, independently of the plot inventory. The watchers exited normally once all tracked jobs were terminal. Temporary SSH interruptions were recovered. RFD4 logged scratch multiprocessing cleanup errors (`Errno 16`) during process finalization; generated CIFs and control/QC results exist and the allocations exited zero, so those messages did not abort inference.

No new GPU job was submitted during this status review. The two-epitope scope and biological validation plan remain unchanged. All **44 local tests**, Ruff and shell syntax checks pass after the future-run limit correction; no measured binding or biological success is inferred.
