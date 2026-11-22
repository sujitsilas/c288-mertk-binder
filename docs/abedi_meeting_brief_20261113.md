# MERTK binder project: meeting with Professor Abedi

November 13, 2026. Scheduler checked at approximately 11:32 a.m. Pacific. This is a discussion brief; no new experiments were submitted during this check.

**Current objective:** Design a MERTK-directed protein that improves efferocytosis. The [functional design plan](efferocytosis_design.md) prioritizes a glycan-compatible Ig1 binding module, followed by a controlled monovalent/bivalent comparison and a conditional apoptotic-cell-targeted bridge. Binding or ligand blockade alone is not success. The lab-reference MERTK models support exploring Ig1, but all six displayed variants overlap the resolved N114 glycan in our structural audit.

**Suggested opening**

“We are designing a MERTK-directed protein to increase macrophage clearance of apoptotic cells. The positive control passed, but we do not yet have a qualified MERTK binding module. Ig1 is our current lead, with N114 glycan access still to solve. I want to compare monovalent and defined bivalent engagement using actual uptake and receptor retention as endpoints, with a PS-targeted bridge as a fallback.”

**Current runs**

| Run | Job | Status | Result |
|---|---:|---|---|
| PD-L1 positive control | 1573755 | Completed, 10m42s | Two accepted sequence variants from one independent backbone; control gate passed. |
| FNIIIb V469/P471/F472 | 1573763 | Stopped at four-hour allocation limit | 30 complete attempts, one interrupted; 280 saved sequence proposals; zero accepted. |
| Ig1 I122/F147/Y148 | 1573761 | Pending, `ReqNodeNotAvail` | Control dependency cleared. Eligible for A100, H200 and 32 GB V100; 12-hour limit, goal of up to 200 attempts. |
| Ig1 native-context QC | 1573762 | Pending on generation | No Ig1 results yet from this new patch. |
| AF3 confirmation | — | Not submitted | No MERTK candidate has been nominated. |

Slurm's start estimate is unstable: one `scontrol` response briefly showed 5:23 p.m. on a V100 node, but the subsequent `squeue --start` returned N/A. There is no reliable start time. About **5.56 of 24 GPU allocation-hours** have been consumed across the project; queued Ig1 time is not consumed time. Two hours remain reserved for AF3. Different GPU models count equally by allocation time for this project cap, not by throughput.

**What the FNIIIb result actually says**

Of 280 saved proposals, 267 failed base AF2 validation. The other 13 reached full scoring, but all derive from one 130-residue backbone. Their mean interface pTM is 0.73–0.79, yet their mean buried unsatisfied polar-atom score is **8–12.5 against a cutoff of 4**. All 37 parent interface positions stayed fixed during ProteinMPNN redesign, so these variants did not explore new interface chemistry.

Separately, their footprints extend across the FNIIIa/FNIIIb junction. All 13 overlap the experimental N294 glycan even after aligning FNIIIb alone: closest heavy-atom separations are 0.31–0.66 Å, below the existing 2.5 Å screen. Glycans were absent from generation, so this did not cause the AF2/Rosetta rejections; it is an additional native-accessibility problem. This comparison uses one resolved glycan conformation and does not establish that every FNIIIb epitope is impossible.

![Rejected FNIIIb design overlapping the N294 glycan](../outputs/analysis/20261113/fniiib_glycan_overlap.png)

Green: rejected binder; orange: N294 glycan; blue/gray: native FNIIIb/FNIIIa; magenta: requested hotspots. This is the most useful structure to bring to the meeting. [Detailed diagnosis and provenance](fniiib_failure_diagnosis.md).

**Four decisions to discuss**

1. **Which Ig1 format is most likely to improve uptake?** The queued I122/F147/Y148 patch tests a different surface from the old pilot. Ask: “Should our first functional comparison use a monovalent parent and a defined bivalent format, and when would a PS-targeted bridge be justified?” Require native access and actual increased clearance. An Ig1 competitor that reduces uptake is an informative failure for this objective.

2. **How should we prevent glycan and domain-junction conflicts during generation?** Ask: “Would you use excluded volumes, ensembles of glycosylated receptor conformations, or a different crop/orientation constraint to keep the whole binder accessible?” The observed FNIIIb problem argues for inspecting the entire binder footprint early, not just a hotspot-centred patch. The structural reference contains modeled carbohydrate at ten of fourteen potential sites and was produced in insect cells; it does not fully represent mammalian glycosylation. [Starbird et al., 2026](https://pmc.ncbi.nlm.nih.gov/articles/PMC13660871/).

3. **Is a failed interface worth rescuing?** Ask: “Once a pose is compatible with the receptor and glycans, should I compare interface redesign against new backbone sampling?” The fixed interface explains why the 13 sequence variants are a narrow test. Interface redesign is possible, but should be a separately recorded comparison with the same evaluation thresholds. It would not repair a binder occupying glycan volume by itself. [BindCraft guidance](https://github.com/martinpacesa/BindCraft/wiki/De-novo-binder-design-with-BindCraft).

4. **How will we demonstrate improved MERTK-dependent efferocytosis?** Ask about primary human macrophages, apoptotic-neutrophil uptake imaging, MERTK dependence controls and a second uptake challenge. Binding/selectivity and GAS6/PROS1 competition qualify the parent; increased internalization with preserved viability and useful receptor retention establishes the desired effect. Agree which experiments are achievable for the course and report an untested candidate as an intermediate deliverable, not functional success.

**Decisions to leave with**

- A primary epitope and one fallback strategy, with a clear reason to switch.
- A practical method to account for receptor glycans and domain context before expanding sampling.
- A small format panel, a donor-paired internalization endpoint, and an achievable stopping rule within the remaining budget.

**Wording to correct when presenting the proposal**

Say “positive-control feasibility demonstrated,” not “published acceptance rate reproduced.” Say “up to 200 attempts within 12 hours,” not “at least 200 guaranteed.” Use the actual Ig1 anchors I122/F147/Y148; Y123 is not in this job. Neither patch is proven glycan-clear. Thirty unsuccessful FNIIIb attempts do not establish that FNIIIb cannot bind, and thirteen correlated variants do not equal thirteen independent binder solutions. Keep the claim at computational design until binding is measured.

**Evidence**

- [Scheduler status](../outputs/hpc/round_20261112_01/status.json) (updated by the monitor; this brief is a dated snapshot).
- [PD-L1 native metrics](../outputs/hpc/round_20261112_01/campaigns/pdl1_control/metrics.json).
- [FNIIIb native metrics](../outputs/hpc/round_20261112_01/campaigns/FNIIIb_F472/metrics.json).
- [Experiment specification, accounting and proposal review](next_experiments.md).
- [Current functional design plan](efferocytosis_design.md).
- [Historical November 12 proposal](repository_archive.md#proposal-drafts).
