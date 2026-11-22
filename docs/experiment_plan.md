# Biological design and experimental validation plan

**Historical program-level plan.** The [November 13 functional design plan](efferocytosis_design.md) defines the current objective and staged panel: improved MERTK-dependent efferocytosis, with a monovalent parent versus defined bivalent Ig1 comparison and a conditional PS-bridge follow-up. The assay concepts below remain useful; the old parallel generator/fusion panel and mouse/tissue expansion are not automatic current requirements.

Draft protocol, updated November 11, 2026 with the mouse-testing requirement. These are proposed experiments and decision criteria, not observed results. The [proposal](repository_archive.md#proposal-drafts) defines the scientific rationale. Lock endpoints, doses, time points, exclusions and analysis after an assay-development pilot and before testing the nominated panel.

## Primary objective and staged scope

The primary endpoint is the donor-paired difference in internalized apoptotic-neutrophil cargo per viable human monocyte-derived macrophage. Record both the percentage of macrophages containing cargo and cargo per macrophage. Increased surface adhesion, reporter brightness or macrophage survival alone does not establish increased efferocytosis.

Stage the work to avoid an unmanageable full factorial screen:

1. Nominate up to four diverse monovalent candidates from each generator, drawn from the same reviewed epitope classes. Test expression, monodispersity and MERTK binding first.
2. Characterize binding hits for AXL/TYRO3 selectivity and compatibility with GAS6 and PROS1. Advance the best two into a defined bivalent format.
3. In parallel, test two Annexin A5–PROS1-module fusion configurations, keeping the natural receptor-binding module intact initially. Confirm module function before considering sequence changes.
4. Test cellular uptake and mechanism in a staged donor assay. A small assay-development pilot establishes dynamic range; a proposed six-donor confirmation study is a feasibility starting point, not a claim of statistical power.
5. Advance successful proteins to resolution assays, followed by a tissue repair model. Tissue regeneration and animal efficacy remain downstream studies.

## Reference and epitope dossier

For each epitope, record the experimental coordinate checksum, human UniProt numbering, author numbering, construct boundaries, omitted residues and glycans, design-chain mapping, proposed contact residues, and all modeling assumptions. The experimentally observed 37NT A93–A484 construct is the shared starting point, not a model of the complete membrane receptor.

| Epitope class | Biological rationale | Required review before generation |
|---|---|---|
| Ligand-compatible extracellular patch | Permit endogenous PS–ligand–MERTK bridging while adding a controllable binding module | Natural-ligand footprint hypothesis, glycan accessibility, MERTK-specific side chains, domain movement and linker exit direction |
| Membrane-proximal FNIIIb patch | Test whether nearby occupancy can reduce receptor shedding | Human stalk/cleavage mapping, accessibility to a soluble binder, membrane collision, ligand effects and perturbation of putative receptor contacts |
| PROS1-derived receptor interface | Preserve natural recognition in an engineered PS bridge | Domain integrity, simultaneous membrane/receptor access, MERTK versus TYRO3 engagement and compatibility with the chosen fusion orientation |

Do not transfer the murine cleavage residue number directly to human MERTK. Resolve human sequence correspondence and, where feasible, identify released fragment termini. A reduction in shedding must be measured; proximity to a predicted cleavage site is not sufficient evidence.

## Human-to-mouse translation requirement

Retain experimental human 37NT as the structural reference, and treat human and
mouse MERTK as positive targets during epitope review and candidate selection.
For planned testing in conventional mice, require demonstrated mouse receptor
engagement and function before animal nomination. Human-only runtime pilots do
not establish mouse applicability. Published anti-MerTK screens found
species-restricted antibodies, including a human candidate without mouse
cross-reactivity ([primary study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10067786/)).

A November 11 alignment of reviewed canonical UniProt sequences
[human Q12866](https://www.uniprot.org/uniprotkb/Q12866/entry) and
[mouse Q60805](https://www.uniprot.org/uniprotkb/Q60805/entry) found:

| Human sequence window | Identity to aligned mouse sequence, counting gaps |
|---|---:|
| Annotated ectodomain, 21–505 | 373/487 columns (76.6%) |
| Observed 37NT segment, 93–484 | 326/393 columns (83.0%) |
| Annotated FNIIIb domain, 386–484 | 74/99 columns (74.7%) |

These are our sequence calculations, not measured binding probabilities. Method:
global alignment of annotated ectodomains with BLOSUM62, gap opening −10 and
extension −0.5; first optimal alignment. Source JSON, entry versions/checksums,
alignment and residue mapping are archived in the ignored
`outputs/species_review_20261111/` directory. The observed 37NT sequence was
verified against human positions 93–484. Review local contact residues, indels,
glycosylation and accessibility instead of choosing an epitope from overall
identity. Verify the planned mouse strain/construct sequence. Any mouse model
derived from the human structure must be labeled as a homology model; it is not
an experimental mouse structure.

Test human and mouse ectodomain binding and native cell-surface recognition
early. Check AXL/TYRO3 selectivity in both species. Before mouse efficacy work,
require improved apoptotic-neutrophil uptake in primary mouse macrophages with
Mertk loss/blockade controls, compatible mouse Gas6/Pros1 signaling, and acceptable
shedding/trafficking. Binding or phosphorylation alone is insufficient. Apply
the same species checks to the Annexin A5–PROS1-derived fusion.

If a useful human candidate lacks mouse activity, a separately characterized
mouse-reactive surrogate can test the mechanism, or an appropriate human-MERTK
knock-in model can be considered. A surrogate's efficacy does not validate the
human candidate itself. A humanized model must express the intended receptor in
the relevant macrophages and retain the ligand/clearance biology being tested.

## Candidate quality and binding

Express proteins in a system compatible with their disulfides and glycosylation requirements. Mammalian expression is a reasonable starting point for the PROS1-derived fusion. Record purity, yield, monomer/oligomer distribution, aggregation and endotoxin burden; contamination could change macrophage behavior. Use matched buffer and tag controls.

Measure concentration-dependent binding with SPR or BLI to human MERTK ectodomain, then confirm recognition of membrane MERTK on cells. Use an assay format that separates monovalent affinity from avidity and nonspecific surface binding. Test AXL and TYRO3 using matched receptor quality, density and detection sensitivity, plus receptor-negative cells and an unrelated extracellular protein control.

A proposed screening target is MERTK binding at or below 100 nM and at least a tenfold monovalent affinity margin over each off-target. These are engineering targets, not literature-derived guarantees or final efficacy thresholds. If off-target binding is undetectable, report the tested range and detection limit instead of assigning an infinite selectivity ratio. Functional selectivity at the active concentration takes precedence over an arbitrary affinity cutoff.

## Natural-ligand competition and apoptotic-cell recognition

Test each candidate in ligand-free, GAS6-supported and PROS1-supported settings, with ligand-alone controls. Use binding-competent, appropriately processed natural ligands and document serum composition; uncontrolled ligand supply can hide competition or create apparent ligand-independent activity. Include sequential-binding measurements in both addition orders and a cell-based functional confirmation.

For Concept A, ask whether candidate binding preserves ligand occupancy and ligand-supported uptake. A loss of ligand binding is a mechanistic concern even if receptor phosphorylation increases. For Concept B, competition at the natural receptor-binding site may occur; measure whether the fusion provides effective PS-dependent bridging and improves net uptake despite that competition.

Assess apoptotic-neutrophil quality with an independent aliquot: apoptosis, membrane integrity, PS exposure and the fraction already undergoing secondary necrosis. Avoid pre-coating cells used for uptake with Annexin V detection reagent, because it can compete for the same surface being tested. Include viable and secondary-necrotic target controls. Confirm PS binding on defined PS-containing versus control membranes; calcium-dependence belongs in a separate binding assay, since chelation in a cellular uptake experiment has broad effects.

For the fusion, compare intact fusion, Annexin A5 alone, PROS1 module alone, unfused modules, full-length PROS1, a validated PS-binding-defective fusion and a validated receptor-binding-defective fusion. Verify that defective variants remain folded. Demonstrate simultaneous engagement rather than two unrelated positive binding assays. Quantify loss of natural-ligand recruitment when PS is heavily occupied.

## Productive receptor organization and signaling

Compare monovalent and bivalent forms of the same binder with matched linker/tag controls. Include equal binding-site concentration comparisons and report avidity separately. Measure receptor proximity or clustering using a suitable imaging/proximity assay, together with MERTK phosphorylation and downstream signaling time courses.

Do not use ligand-free soluble ectodomain aggregation as proof of productive cell-surface clustering. Rule out nonspecific protein aggregation, Fc receptor effects and reporter overexpression artifacts. Test signaling in the presence and absence of apoptotic cells to determine whether activation depends on the intended contact. Use clearance and second-challenge function to decide whether an observed cluster is useful.

## Shedding, internalization, recycling and degradation

Compare resting and inflammatory conditions that produce measurable MERTK shedding while maintaining macrophage viability. Sample matched time courses for:

- Surface MERTK with a noncompeting detection reagent.
- Soluble MERTK in supernatant with spike-recovery controls.
- Total MERTK and cell-associated cleavage products by an orthogonal assay.
- Receptor internalization, recycling after washout, and degradation/lysosomal routing.
- Cell number, viability and receptor transcript abundance when needed to distinguish expression changes.

Use an ADAM17 perturbation as a mechanistic control and a cleavage-resistant receptor comparator only after the construct and species mapping are justified. Broad ADAM17 effects on other substrates make an inhibitor response insufficient on its own. Detection-epitope masking, sequestration of soluble MERTK, altered transcription or receptor internalization can each mimic reduced shedding.

Expose cells to a second apoptotic-neutrophil challenge after initial uptake. Productive receptor internalization may be normal; persistent surface depletion and loss of repeated clearance are the failure modes. Test soluble MERTK as a potential candidate sink. Claim direct steric protection only if the binding/cleavage evidence supports it; otherwise describe the measured effect as preservation of surface receptor with mechanism unresolved.

## Efferocytosis assay and causal controls

Use primary human macrophages from independent donors and standardized apoptotic-neutrophil preparations. Pair treatments within donor and balance target-cell batches across conditions. Include vehicle, nonbinding protein, each natural ligand, candidate alone and candidate plus each ligand. Use blinded imaging or flow analysis with extracellular fluorescence quenching/stripping and orthogonal image confirmation to separate engulfed cells from attached cells. Normalize to viable macrophages and report target-cell availability.

Establish MERTK dependence with receptor loss or blockade and an orthogonal perturbation. A rescue experiment with functional MERTK is particularly useful; a kinase-impaired comparator can test signaling dependence. A small-molecule inhibitor alone is insufficient because off-target kinase activity could affect uptake. Assess AXL/TYRO3 activation at the functional dose, and distinguish residual endogenous uptake from the candidate-specific benefit.

Measure repeated uptake, digestion/processing where feasible, and residual apoptotic-cell burden. A candidate that only increases adhesion, damages viable targets, or suppresses uptake with endogenous ligand present does not meet the goal. A receptor-protective candidate may be beneficial without being independently agonistic; label its mechanism accordingly.

## Resolution and tissue extension

In a defined inflammatory macrophage–neutrophil model, follow extracellular apoptotic cells, secondary necrosis and inflammatory mediator trajectories after stimulus removal. Include candidate effects on neutrophil viability/apoptosis so accelerated disappearance is not mistaken for improved clearance. If pro-resolving lipid mediators are measured, use validated analytical methods and appropriate standards rather than an isolated surrogate marker.

Advance to a macrophage-containing human skin wound co-culture or explant. Compare vehicle, natural ligand, active candidate, nonbinding/module controls and MERTK-disrupted macrophages. Include tissue plus candidate without macrophages, and control carryover in conditioned-medium experiments. Readouts should include re-epithelialization, barrier function, tissue viability, architecture and matrix/scar formation. Separate effects on cell proliferation, contraction and migration from restoration of tissue.

A future sterile-inflammation model should measure neutrophil numbers over time, distinguishing a lower recruitment peak from faster clearance after the peak. A resolution interval can be estimated from peak to half-peak only with adequate temporal sampling. A separate wound model evaluates tissue repair; neither study automatically demonstrates regeneration. Verify species cross-reactivity before animal work and use the applicable human-sample/animal protocols.

## Analysis and advancement rules

Use donors as the biological replication unit; technical wells are nested measurements. Analyze paired treatment effects with donor blocking or a mixed-effects model, report effect sizes and uncertainty, and retain all prespecified exclusions and failures. Include neutrophil batch as a blocking factor where the design permits. Lock the primary contrast and endpoint time/window before the confirmation cohort; correct multiple candidate comparisons and designate mechanism readouts as secondary.

Choose confirmation sample size from pilot donor variability and a prespecified biologically meaningful uptake improvement. A proposed 20% relative improvement over the matched ligand condition can be used for planning, but it is a hypothesis-dependent threshold to finalize before the confirmatory study. If the baseline is near zero, use an absolute difference. An inconclusive interval is not evidence of equivalence.

Advance candidates only when the uptake benefit is reproducible, linked to MERTK, and compatible with viable-cell discrimination and sustained receptor function. Keep binding, agonism, shedding protection, efferocytosis, resolution and tissue effects as separate labels in the candidate dossier. Report negative results as informative tests of epitope, valency or bridge design.

## Course milestones

| Milestone | Deliverable | Decision |
|---|---|---|
| Reference and epitope review | Glycan-aware 37NT dossier, human/mouse mapping, negative-target panel | Proceed only with an accessible, biologically justified patch |
| Computational pilot | Common-budget generation logs and structural review | Nominate diverse proteins; do not lower filters to fill quotas |
| Biochemical pilot | Expression, monodispersity, MERTK/off-target binding and competition | Advance binding hits with interpretable mechanism |
| Cellular validation | Donor-paired uptake, trafficking and shedding evidence | Advance reproducible functional hits |
| Downstream extension | Temporal resolution and tissue-function measurements | Make only the claims directly supported |
