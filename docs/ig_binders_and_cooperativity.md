# Ig-domain binders and the two-epitope experiment

> Historical scope: the revised proposal now accepts either a non-blocking probe or an antagonist and allows up to six nominees. The [current run plan](next_experiments.md) supersedes the old epitope choices, four-candidate limit and ligand-compatibility pass/fail rule. Historical experiments remain unchanged.

November 12, 2026. This user-requested extension adds one exploratory Ig1 surface to the existing FNIIIa Patch A campaign. The class question remains specific MERTK binding with natural-ligand compatibility. At most four experimental nominees are selected across both epitopes and both generators; the shared 24 allocated A100-hour cap remains unchanged. The existing submission PDF predates this extension.

## Why add the Ig region?

Ig domains are targetable protein surfaces. They require care here because natural ligand recognition involves this region. An Ig binder could coexist with GAS6/PROS1, obstruct them, or alter their binding indirectly. None of those outcomes follows from its domain label alone. Functional mouse Ig mutants affect ligand-dependent retinal phagocytosis, but are not an atomic footprint for human MERTK in macrophages ([Parinot et al., 2024](https://doi.org/10.3390/ijms25126630)). The experimental receptor structure supplies an unliganded reference, not a receptor–natural-ligand complex ([Starbird et al., 2026](https://doi.org/10.1042/BCJ20260584)).

Retain Patch A's receptor-protection motivation as a hypothesis. Binding I365/W367/A369 does not establish inhibition of ADAM17 cleavage. The missing membrane-proximal stalk is not modeled by 37NT. Measure soluble ectodomain release together with surface and total receptor before claiming protection; receptor internalization can also reduce surface signal. Neither an Ig binder nor a Patch A binder can reconnect an already shed ectodomain to the membrane kinase.

## Provisional Ig1 side surface

Use **K178/E183/V185**, with the complete observed Ig pair and intervening linker, **37NT A93–A285**, retained for design. The crop has 193 native residues. Both generators retain native receptor sequence and use an 80-residue monovalent binder. Restore the full A93–A484 protein and its observed glycans for candidate screening. These are design anchors, not experimentally established energetic hotspots.

| Human anchor | Existing solvent-accessible area, Å² | Nearest observed glycan heavy atom, Å | Aligned mouse residue |
|---|---:|---:|---|
| K178 | 129.1 | 17.9 | K173 |
| E183 | 116.7 | 15.5 | E178 |
| V185 | 71.2 | 12.1 | V180 |

Values derive from the existing 37NT PyMOL solvent-area inventory and original mmCIF coordinates. Anchor C-alpha separations are 4.2–8.8 Å. Conservation does not establish mouse binding; neighboring sequence, glycoforms and receptor conformation remain relevant. Observed insect-cell glycans underrepresent mammalian glycan volume and heterogeneity.

For a cautious ligand comparison, align the Ig sequences of experimental GAS6–AXL [2C5D](https://www.rcsb.org/structure/2C5D), then rigid-fit AXL chain C to the MERTK Ig1 C-alpha correspondences. Use **both GAS6 chains A and B**, because the crystallographic assembly contains more than one receptor–ligand interface. The 90-correspondence Ig1 fit has 2.58 Å C-alpha RMSD. Minimum anchor-heavy-atom distances to either projected GAS6 are approximately 11.6, 15.8 and 16.6 Å for K178, E183 and V185, respectively. These distances favor a side-surface trial over deliberately targeting the inferred ligand face; they do not guarantee room for an entire binder or predict PROS1 compatibility. A predicted binder must undergo a full-footprint steric review.

The hypothesis is weaker than a directly solved MERTK–ligand footprint. Human T145/F147 align to mouse T140/F142 functional mutants; human K268/K274 align to mouse K263/K269. Mouse G122 is an alignment insertion without a direct human residue counterpart. Avoid treating mutant positions as a complete, transferable contact map. The chosen anchors avoid deliberately centering those positions. AXL/TYRO3 binding and natural-ligand compatibility remain experimental selection criteria, not assumed properties.

In the same sequence alignment, human MERTK K178/E183/V185 correspond to AXL F113/T118/V120. Two differing anchor identities make selectivity plausible to explore, but do not establish it; a binder can contact a wider, conserved surface. TYRO3 and AXL must both remain negative-target assays.

## Co-binding first, cooperativity second

There are two distinct compatibility questions: can a Patch A binder and an Ig binder occupy the same receptor, and can each binder coexist with GAS6 or PROS1? Test natural ligands separately. Sequential SPR/BLI in both addition orders can identify blocking and potential co-occupancy. Include receptor alone, each component alone, unrelated binder and reference-surface controls; confirm the best ternary complex by an orthogonal solution method, such as SEC-MALS if the complex is stable enough. A mass increase alone does not prove a specific 1:1:1 complex.

**Cooperativity means an affinity change**, not simply simultaneous occupancy. For two monovalent binders A and I and a monomeric receptor R, define

`alpha = Kd(I binding R) / Kd(I binding R·A)`.

Alpha greater than one indicates positive coupling; approximately one indicates independent binding; below one indicates negative coupling. The reciprocal experiment, A binding R versus R·I, must be consistent with the same equilibrium coupling under the model. Fit concentration-dependent equilibrium data to the full ternary binding model, accounting for free concentrations and depletion, rather than calculating alpha from a single sequential trace. Report uncertainty; an interval spanning one is inconclusive. A two-ligand thermodynamic treatment is described in [Velazquez-Campoy et al., 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC1544317/).

For a class project, measure individual affinity and ligand compatibility first. Only if both epitopes yield confirmed monomeric binders, test one best cross-epitope pair using full titrations with and without near-saturating partner, reciprocal addition order, and at least three independent assay runs. Apply the same logic to binder–GAS6 or binder–PROS1 coupling only after nonblocking binding is established. Do not require positive cooperativity for a successful binding-probe project.

Use monovalent reagents and a verified monomeric receptor, with low-density immobilization and checks for aggregation/mass transport. Bivalent constructs, Fc fusions and dense receptor surfaces introduce **avidity**, which can imitate affinity enhancement. Hill slopes and stronger combined signaling do not, by themselves, demonstrate thermodynamic cooperativity. AF3 can suggest steric compatibility but cannot measure alpha or establish productive receptor clustering; the current repository scorer handles two chains and does not implement this ternary analysis.

## Keep function and shedding claims separate

If a pair passes binding and compatibility tests, compare neither binder, Patch A alone, Ig alone and the mixture at matched individual concentrations. Include natural ligand alone and the corresponding ligand-supported conditions, a MERTK pathway control, viability, surface receptor, phosphorylation and true internalized apoptotic-cell uptake. Exclude surface-adherent material from the uptake readout. A mixture exceeding either single binder is not sufficient to claim synergy: specify an additivity model and uncertainty in advance.

Shedding protection is an optional mechanistic follow-up with matched inflammation/protease conditions, soluble MERTK, surface receptor and total receptor. Keep fusion design, PS-bridging engineering, forced clustering and mouse efficacy outside this first two-epitope pilot. Independent ligand-compatible binding at both surfaces is already a useful outcome.

## Run specification

Fresh campaigns preserve the failed Patch A input-check run. Each new epitope receives at most two BindCraft trajectories/four A100-hours and two RFD4 replicates/two A100-hours, with a separate one-replicate upstream PD-L1 control/one A100-hour gate. The two campaigns reserve at most 14 A100-hours, leaving ten within the existing cap for assessment or explicitly accounted follow-up. No large screening expansion is implied. Failure of the RFD4 control blocks that campaign's RFD4 MERTK job. Input validation runs first; generation is submitted only after its success is observed.

## Local structural review

Reproduce the anchor measurements with `.venv/bin/python scripts/review_ig_epitope.py`. Render with the local Python providing `pymol2`, for example `/opt/homebrew/opt/python@3.14/bin/python3.14 scripts/render_ig_epitope.py`. The default output directory is `outputs/epitope_review_20261112/`: `review.json` records input hashes, sequence-alignment settings, residue correspondences and the rigid transform; `MERTK_Ig1_review.pse` contains overview and Ig1 close-up scenes. Both exported PNGs were visually inspected. The native receptor and observed glycans are shown without an invented natural-ligand complex. A copy is saved in `~/Downloads/MERTK_Ig1_Review_2026-11-12/`.
