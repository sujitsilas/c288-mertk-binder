# Abedi meeting follow-up: glycan avoidance and MERTK cleavage protection

**November 17 update:** A separately authorized cleavage-site stalk pilot is now submitted. The November 13 scope notes below are history; the old FNIIIb F472 campaign remains retired. [Current protector design](cleavage_protector_20261117.md).

**Superseded scope, November 13 evening:** The user subsequently retired FNIIIb and restricted every new binder epitope to Ig1/Ig2. The direct membrane-proximal/stalk design proposed below is historical and will not be launched. The active protector arm screens monovalent Ig1/Ig2 parents for distal protection; full receptor/stalk context remains useful for interpretation. See the [current design plan](efferocytosis_design.md) and [Ig-only target inputs](../outputs/design/ig_only_20261113/manifest.json).

November 13, 2026. The project objective remains improved MERTK-dependent efferocytosis. This note records the professor's feedback, the current implementation and a second design mechanism.

**Live run snapshot, 6:11 p.m. Pacific**

| Experiment | Job | State and observed result |
|---|---:|---|
| Ig1 I122/F147/Y148 | 1573761 | Running on g9, Tesla V100-PCIE-32GB. Started 1:25:50 p.m.; 12-hour limit ends November 14 at 1:25:50 a.m. unless it finishes earlier. |
| Ig1 native-context report | 1573762 | Pending on generation. A separate interim local check was run now. |
| FNIIIb V469/P471/F472 | 1573763 | Four-hour timeout; 30 complete attempts, zero accepted variants. |
| PD-L1 positive control | 1573755 | Completed; two accepted variants from one backbone; control gate passed. |
| AF3 / cleavage-protector generation | — | Not submitted. |

The synchronized Ig1 snapshot contains **nine starts, eight complete attempts and one in progress**. Seven trajectories reached structural scoring. Two 123-residue variants, `Ig1_Q124_l123_s483444_mpnn1` and `mpnn3`, pass BindCraft filters and belong to one backbone. They have mean interface pTM 0.76/0.77, normalized interface PAE 0.32/0.32 and buried-unsatisfied-polar scores 1.5/3.5. These are computational filter passes, not measured binding or nominated candidates. Counts can increase while the run continues.

Both fail native-context screening. After fitting only Ig1, their target C-alpha RMSDs are 0.58 and 0.56 Å, but closest binder–N114-glycan heavy atoms are 0.65 and 0.65 Å. Thus the glycan conflict persists after avoiding a whole-Ig-pair alignment. The generated target is renumbered A1–A193; it was mapped by exact sequence to native A93–A285 before selecting Ig1. This prevents confusing generator numbering with human residue numbering. No candidate is nominated for AF3.

Evidence: [live metrics](../outputs/analysis/20261113/ig1_live_metrics.json), [native-context screen](../outputs/analysis/20261113/ig1_live_native_context.json), [domain alignment checks](../outputs/analysis/20261113/ig1_live_alignment_sensitivity.json). The generic metrics summarizer labels an unfinished attempt as censored; at this snapshot it is still in progress, not a terminated attempt.

**Are glycosylated sites already protected from binder overlap?**

Partly. Here, “protect” means reserve glycan-occupied space so that the designed binder does not obstruct it. It does not mean that the binder must cover a glycan attachment site.

| Step | Current implementation |
|---|---|
| Site prioritization | All 14 potential N-linked sequons contribute to a heuristic exclusion-radius score; observed glycan distances are recorded. This penalizes risky patches, but is not a hard constraint. |
| Running BindCraft generation | Protein-only target; hotspots are supplied. No glycan atoms, coldspot repulsion or glycan excluded-volume loss is used. The frozen manifest explicitly records `glycans_modeled=false`. |
| Native-context screening | Accepted poses are checked against intact protein and resolved glycans; the observed-glycan distance gate is 2.5 Å. Unresolved sites and mammalian glycoforms need further review. |
| Interim diagnosis | All six lab-reference Ig1 variants and both current accepted Ig1 variants overlap N114 carbohydrate in the examined native overlays. |

The next generation experiment should therefore add **during-design avoidance**, followed by the existing whole-pose checks. Defining a glycosylated Asn as a coldspot is useful but does not by itself exclude the carbohydrate's entire volume. Preserve glycan-bearing receptor sequence, consider the surrounding surface and plausible glycan conformations, and reject any pose that requires an artificial crop edge or removal of native carbohydrate.

The current pinned BindCraft engine has no configured coldspot option. BindCraft2 documents target `coldspots` and a repulsion objective, making it a candidate implementation route; a new pinned environment and control validation would be required before using it. It is not the engine running this job. Explicit carbohydrate-volume avoidance and post-design screening remain necessary regardless of residue masks. [BindCraft2 settings](https://github.com/PacesaLab/BindCraft2/blob/main/docs/source/reference.md).

**Second mechanism: a MERTK-specific cleavage protector**

Design a **monovalent MERTK-binding protein that reduces ectodomain shedding**, preserving intact surface receptor for endogenous ligand-supported clearance. This mechanism can improve function without forced receptor crosslinking. It is an explicit second design arm; the existing Ig1 and FNIIIb jobs do not test it.

There is direct precedent: Duprez et al. reported MERTK-selective nanobodies that inhibited shedding while preserving GAS6 binding. Their bivalent format reduced surface receptor. The study leaves direct steric shielding versus allostery unresolved and does not establish increased efferocytosis by our proposed molecule. [Duprez et al., 2026](https://pmc.ncbi.nlm.nih.gov/articles/PMC13378718/). Cleavage-resistant Mertk improved efferocytosis in an atherosclerosis model, supporting the biological hypothesis through a separate genetic intervention. [Cai et al., 2017](https://jci.org/articles/view/90520).

For a direct shielding design, prioritize the membrane-proximal FNIIIb/stalk region rather than assuming an Ig1 binder will protect cleavage. The 2026 paper reports the human P490–S491 cleavage bond. Those residues are beyond the last observed 37NT residue, G484, so our present structural crops cannot represent the bond. The frozen UniProt reference supplies extracellular residues 485–505, `WVDYAPSSTPAPGNADPVLII`, but their conformations must be modeled with uncertainty and membrane context. The earlier F472 footprint is not established as protective. [Target-preparation specification](../outputs/analysis/20261113/cleavage_protector_target_spec.json).

Next preparation steps are to build and review a membrane-proximal conformational ensemble, preserve nearby glycan access, choose a surface that could obstruct protease approach without obstructing ligand function, and then specify a separate bounded campaign. A peptide-only binder to the cleavage sequence would also require verification on full-length membrane MERTK.

The protection experiment needs concordant **less soluble MERTK, preserved intact surface receptor, preserved ligand responses, and improved internalized apoptotic-cell uptake** under a cleavage-promoting condition. Include total receptor and orthogonal detection to exclude epitope masking or soluble-MERTK assay interference. Test GAS6 and PROS1 separately; assess whether the binder inhibits ADAM17 generally or protects MERTK specifically. A known shedding inhibitor or cleavage-resistant receptor is a mechanistic comparator, not evidence of binder specificity. Normal basal uptake must remain intact.

The receptor-preserving route is a credible alternative to the bivalent agonist route. Both must pass the same efferocytosis endpoint; neither binding nor reduced shedding alone is sufficient. No new jobs or changes to the frozen running experiment were made during this review.
