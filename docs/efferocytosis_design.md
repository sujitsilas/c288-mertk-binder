# Two MERTK binder-design aims for improved efferocytosis

The [November 18 two-page proposal](revised_project_direction.md) uses pooled protease transcript counts from the mouse cohort to prioritize tests across inflammatory settings. It also includes the Hoffman2 test of the Ig1 binders as direct shedding masks. The November 17 design history follows below.

The user's clarified project has **two separate deliverables**:

| Aim | Molecule | Evidence required |
|---|---|---|
| 1. Ig-domain binder | A selective, glycan-compatible binder to MERTK Ig1 or Ig2 | Native binding and accessibility first; receptor-pairing and apoptotic-cell engagement formats then test whether it improves uptake. |
| 2. Cleavage protector | A MERTK-binding molecule that reduces receptor shedding | Reduced released ectodomain together with retained intact surface receptor and preserved ligand-supported efferocytosis. |

The two products may differ. A shared Ig-domain parent can enter both tracks, but a result in one does not establish the other. Dimerization and PS targeting are functional extensions of Aim 1, not additional project aims or prerequisites for the monovalent protection screen. [Machine-readable aims](../configs/project_aims.json) · [Protector design specification](../configs/mertk_cleavage_protector.json).

Current direction, November 17, 2026. Aim 1 remains restricted to **Ig1 or Ig2**. The user has now requested independent triage and execution of Aim 2; that pilot targets the **P490–S491 cleavage-site stalk**, with the old FNIIIb F472 epitope still retired. The objective is macrophage internalization of apoptotic cells, with receptor pairing, cargo recruitment and shedding protection treated as mechanisms to evaluate. No current candidate has demonstrated binding, dimerization or improved efferocytosis experimentally.

**Operational scope**

The November 14 Ig suite ended with two new Ig1 variants failing native protein geometry and no accepted Ig2 variant. Both comparator AF3 jobs completed. Its protector step was only CPU triage. The November 17 suite now includes **independent stalk-targeted generation**, GPU integration preflight, conformer/context screening and conditional AF3. [Protector triage and submitted jobs](../docs/cleavage_protector_20261117.md).

**Binding modules and glycan constraints**

The previous Ig1 run reached its 12-hour limit with 23 completed attempts and one interrupted attempt. Eight variants from four backbones passed BindCraft filters, but none passed native glycan checks. The new experimental loss penalizes overlap during generation, subject to a successful GPU forward/backward preflight; it does not establish glycan-compatible binding in advance.

Gokul's supplied reference contains six displayed Ig1 variants from two backbone families. It supports pursuing Ig1, but the audited poses all overlap N114 carbohydrate. The final baseline contains eight filter-passing variants from four Ig1 backbone families; all overlap N114 under local Ig1 alignment. Neither set is a qualified parent for scaffold assembly. The reference supplies no equivalent evidence establishing Ig2 performance.

| Prepared site | Preferred anchors | Rationale and unresolved issues |
|---|---|---|
| Ig1 E152-centred patch | I122, Y148 | A narrower approach near the reference Ig1 footprint. The heuristic patch has zero members within 15 Å of a sequon CA; that does not establish whole-binder clearance. F147 is removed as a preferred anchor, but contacts to it and ligand competition remain possible. |
| Ig2 L270-centred patch | L270 | A second Ig-domain hypothesis. L196 is excluded by the conservative N170 mask. Higher-ranked I283/S285 patches extend into FNIIIa and are outside the active scope. The remaining patch is still ligand-proximal and glycan-risky. |

Both use native A93–A285 Ig-pair coordinates and preserve residue numbering. Whole-receptor context remains mandatory when evaluating a pose. Retiring an epitope does not justify removing that domain from steric checks. The bundle includes all 14 potential N-linked sequons, 420 resolved carbohydrate heavy atoms and explicit generator-to-native residue mapping. The provisional residue mask flags CA distances below 15 Å from sequons or heavy-atom distances below 6 Å from observed sugars. These are conservative review heuristics, not validated mammalian glycan envelopes. The later experimental generation loss is described in the current queue document.

[Prepared target bundle](../outputs/design/ig_only_20261113/manifest.json) · [Glycan constraints and site audits](../outputs/design/ig_only_20261113/glycan_exclusion_review.json) · [Architecture specification](../configs/ig_scaffold_design.json).

**A concrete scaffold comparison**

The lead architecture hypothesis is **two copies of one qualified MERTK-binding module plus one apoptotic-cell recognition module**. Start with an Fc-free scaffold to keep the uptake mechanism attributable to MERTK. Receptor dimerization is the intended action; making the scaffold itself dimeric does not establish productive receptor pairing.

| Construct | MERTK-binding arms | PS-binding modules | Question |
|---|---:|---:|---|
| M1 | 1 | 0 | Does the parent bind native MERTK, preserve ligand function and alter uptake? |
| M2 | 2 | 0 | Does engaging two receptor molecules improve uptake? |
| M1–PS | 1 | 1 | Does cargo recruitment improve uptake without an added receptor-binding arm? |
| M2–PS | 2 | 1 | Does combined receptor pairing and apoptotic-cell recruitment outperform the matched simpler formats? |
| P1 | 1 | 0 | Does a selected Ig1/Ig2 parent preserve surface MERTK during shedding? |

Use the same parent and cargo module across comparisons. Choose linkers and domain order from qualified bound structures, terminal exit directions and feasible membrane geometry. No exact linker or assembled sequence is justified yet. A tandem Ig1/Ig2 binder could occupy both domains on one receptor; it cannot be counted as a dimerizer without evidence for two distinct MERTK molecules.

Human annexin A5 is an available project starting point for phosphatidylserine (PS) recognition, with module selection still open to evidence. PS masking can counteract uptake, so demonstrate simultaneous receptor/cargo binding and test the PS module alone, an unfused mixture and folded binding-defective controls. Natural GAS6 and PROS1 remain functional benchmarks: ligand binding, PS engagement and signaling are separable properties. [Lew et al., 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4206827/).

MerTK-directed bispecifics have produced targeted phagocytosis in other cargo settings. That supports testing receptor–cargo coupling, while our proposed apoptotic-cell construct remains unvalidated. [Kedage et al., 2020](https://doi.org/10.1080/19420862.2019.1685832).

**Dedicated cleavage-site protection**

**First establish the benefit of cleavage resistance.** Compare wild-type and experimentally confirmed cleavage-resistant human MERTK under matched resting and inflammatory conditions. Match baseline receptor surface abundance and confirm GAS6 responsiveness. Read out released ectodomain, intact surface receptor, and apoptotic-cell uptake per viable macrophage. A reduction in shedding with improved uptake supports testing an unmasked blocker. If shedding falls and uptake remains unchanged, check whether shedding limits uptake in this assay before committing to further design. Mouse studies provide a precedent for the causal test. [Thorp et al., 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3190938/) · [Cai et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4988577/).

The exploratory screen generated monovalent poses around Y488/P490/S491 using three modeled stalk conformations. The 3.41 Å experimental structure ends at G484. AFDB v6 confidence averages 43.4 across residues 485–505. This segment may be flexible and proximal O-glycosylation could influence access to the cut. Its glycoform occupancy remains unknown. The single-sequence Boltz-2 run failed to recover a high-confidence GAS6 control in the expected Ig region. The modeled screen identifies poorly fitting conformations and provides limited evidence about binding. A future control run would test natural-ligand MSAs, including pairing where supported, and judge recovery of the expected Ig-region geometry alongside the scores. [37NT](https://www.rcsb.org/structure/37NT) · [TAM ectodomain structures](https://pmc.ncbi.nlm.nih.gov/articles/PMC13660871/) · [O-glycosylation and shedding](https://pmc.ncbi.nlm.nih.gov/articles/PMC4664366/).

Published MERTK nanobodies inhibited shedding while preserving GAS6 binding. Their bivalent format increased ERK signaling and reduced surface receptor. If the cleavage-resistance test supports intervention, prioritize discovery or screening on native glycosylated MERTK and use the stalk peptide to bias selection when useful. Check hits on cell-surface receptor and test whether they distinguish intact receptor from shed ectodomain. This can capture protective binding modes beyond the uncertain modeled stalk poses. [Duprez et al., 2026](https://pmc.ncbi.nlm.nih.gov/articles/PMC13378718/).

Require concordant reduction in soluble MERTK, retention of intact surface receptor, preserved total receptor and ligand-supported uptake. Use independent detection epitopes and controls for soluble-protein assay interference, receptor internalization, viability and broad ADAM17 inhibition. Test GAS6 and PROS1 separately. Combine a proven protector with a functional bridge only after measuring the individual mechanisms and recheck the complete fusion.

**Unmasked blocker and shedding reporter**

Use the same ELISA and intact surface MERTK measurements for any unmasked blocker that advances. Compare stimulated vehicle, blocker, and unrelated protein alongside unstimulated baseline and cleavage-resistant receptor. The cap and GPQG↓IWGQ linker are an optional follow-up for inflammation-dependent activation. Cleavage of that tether releases the cap and exposes the blocker. [MerTK nanobody shedding assay](https://pmc.ncbi.nlm.nih.gov/articles/PMC13378718/) · [Protease-removable miniprotein masks](https://pmc.ncbi.nlm.nih.gov/articles/PMC12703739/).

For tagged confirmation, pilot extracellular HiBiT on the ectodomain side of P490–S491, away from the stalk and ligand-binding surface. Measure HiBiT luminescence in cleared conditioned medium and on washed intact cells in paired wells. Confirm the expected released fragment, track viability and validate GAS6 binding, surface expression and shedding against untagged MERTK. A C-terminal fluorophore can help check expression. Use the extracellular tag or a surface antibody to quantify surface receptor because the C-terminal fragment undergoes subsequent processing. [HiBiT reporter](https://pubmed.ncbi.nlm.nih.gov/28892606/) · [MERTK fragment study, Figure 4](https://www.mdpi.com/1422-0067/25/8/4404).

**Evidence needed to advance**

First qualify a parent through whole-pose glycan/accessibility review, independent structural assessment, expression, aggregation/state assessment, native glycosylated MERTK binding, AXL/TYRO3 selectivity and ligand-competition measurements. Then establish actual receptor stoichiometry and cell-surface proximity for M2/M2–PS. A high-confidence complex prediction cannot establish productive signaling.

Retain E152 and the 114-residue Q124 design as distinct Ig1 starting backbones. Q124 recovered 23 of its 25 source contacts with 5.62 Å contact PAE and 0.80 ipTM. Its experimental overlay flags protein and glycan interpenetration, so assess placement and glycan compatibility before nominating it. Refolding AF2-selected designs with Boltz-2 measures recovery with a related predictor. Direct binding to glycosylated MERTK is the functional decision point.

The primary endpoint is a donor-paired increase in **internalized apoptotic-cell cargo per viable human macrophage**, with appropriate MERTK dependence. Distinguish internalization from adhesion, record ligand/serum conditions, compare matched binding-site exposure, and assess surface/total receptor and repeated uptake. Include viable-cell uptake and matched nonbinding-scaffold controls. Choose confirmatory donor numbers and effect margins from assay variability rather than treating technical wells as independent samples.

An increase in phosphorylation or dual-surface binding is an intermediate result. A construct that depletes receptors and reduces later clearance fails the intended sustained function. The desired response is effective efferocytosis with useful receptor retention; dimerization and PS binding are hypotheses toward that outcome.

**Current deliverable and remaining work**

Ig1/Ig2 generation and comparator AF3 have run. A separate cleavage-protector design campaign is submitted with glycans and a modeled membrane exclusion term. No native-binding parent, protective sequence, dimerizing scaffold or improved efferocytosis is established; structural scores are hypotheses for experimental testing.
