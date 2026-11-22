# MERTK, GAS6 and PROS1: local PyMOL review

Prepared November 11, 2026. The local review bundle is `outputs/pymol/mertk_review_20261111/`, with a copy in `~/Downloads/MERTK_PyMOL_Review_2026-11-11/`. Open **Open_in_PyMOL.command** from the bundle, or run `pymol outputs/pymol/mertk_review_20261111/MERTK_ligands_hotspots.pse` from the repository. All coordinates are embedded in the session and also available under `structures/`. No server connection is needed to view the session.

The [November 12 first-principles review](first_principles_bridge.md) adds a separate 12-scene session with PS chemistry, experimental calcium/Gla/lipid and Annexin/headgroup complexes, and the distinction between protecting intact MERTK and recovery after shedding. It does not replace this bundle's original failure overlays.

## Start with the saved scenes

Click a scene button in PyMOL, or enter `scene 03_patch_A, recall` in its command line.

| Scene | What to inspect |
|---|---|
| `01_MERTK_domains` | Experimental human 37NT. Ig1 cyan, Ig2 blue, FNIIIa teal, FNIIIb salmon; observed glycans gold. |
| `02_three_molecules` | Left: MERTK; center: experimental GAS6 LG1/LG2; right: predicted PROS1. Molecules are separated and independently oriented for display. |
| `03_patch_A` | Proposed FNIIIa anchors I365, W367, A369 in magenta. |
| `04_patch_B` | Exploratory FNIIIb anchors R460, P471, F472 in magenta. |
| `05_mouse_changes` | Human residues that differ from mouse shown in magenta; remaining protein gray. |
| `06_AXL_GAS6` | Experimental GAS6–AXL complex 2C5D; an analogy for ligand recognition, not a MERTK complex. |
| `07_PROS1_domains` | Predicted PROS1 monomer. Gla green, EGF blue, LG1 yellow-orange, LG2 orange. |
| `08_PROS1_confidence` | PROS1 pLDDT, red low to blue high. |
| `09_GAS6_full` | Predicted full GAS6, same domain colors. |
| `10_PROS1_EGF_NMR` | Experimental PROS1 EGF3–4 fragment 1Z6C, state 1 of 20 NMR models. LG1/LG2 are absent. |
| `11_RFD4_rep0`, `11_RFD4_rep1` | Experimental MERTK gray, generated MERTK colored, designed chain green. Inspect distortion after receptor alignment. |

For the first patch, `show sticks, A_FNIIIa` and `zoom A_FNIIIa, 15` are useful. `show surface, MERTK_37NT` displays receptor surface. `set transparency, 0.5, MERTK_37NT` makes it transparent. `show spheres, unmodeled_glycan_sites and name ND2` marks potential glycan sites whose sugars were unresolved.

## Proposed patches, not established energetic hotspots

Patch A is the better first surface to review: it lies on FNIIIa, away from the Ig ligand-binding module. I365/W367/A369 are solvent exposed, conserved in canonical mouse MERTK, and differ at the sequence-aligned positions of human AXL and TYRO3. These are approximate paralogue sequence correspondences, not measurements of selectivity. Nearby N316/N354 glycans may obstruct a larger binder footprint. Compatibility with GAS6/PROS1 still needs testing.

Patch B explores an alternative FNIIIb surface. Its anchors are conserved in mouse, but several are also similar or identical in AXL/TYRO3. It offers weaker sequence-based discrimination. FNIIIb can participate in receptor contacts; occupancy could alter clustering. Its location does not establish protection from stalk cleavage.

`candidate_patches.csv` contains the human and mouse numbering, measured solvent exposure, distance to observed glycan atoms, and approximate paralogue correspondence for each anchor. `review.json` records source hashes, methods, limitations and failure diagnostics. Neither patch has been submitted to a generator.

## What the structures establish

- **MERTK 37NT:** experimental human ectodomain, observed residues 93–484. The complete juxtamembrane stalk, membrane and kinase are absent. Insect-cell glycans are incomplete and differ from mammalian glycoforms. Fourteen N-X-S/T sequons are marked as potential sites, not assumed occupied.
- **GAS6 1H30 / AXL–GAS6 2C5D:** experimental human receptor-binding fragments. The 2C5D receptor is AXL. It cannot establish a MERTK binding pose.
- **PROS1 and full GAS6:** downloaded AlphaFold DB v6 predicted monomers, not predictions run by this project and not experimentally determined MERTK complexes. Signal/propeptide segments are hidden in domain views; interdomain placement remains uncertain. Native gamma-carboxylation, complete glycans and the PS membrane are not represented.
- **PROS1 1Z6C:** experimental EGF-domain fragment, which does not cover the receptor-binding LG module.

## Pipeline diagnosis

BindCraft's original Conda setup failed and canceled its generation dependency. Recovery required compatible CUDA wheels and an NVCC 12.6.77 pin for JAX 0.4.35. Recovery installation 1459732 succeeded, and A100 runtime check 1462338 passed on November 12. This verifies GPU imports, JIT matrix/convolution operations and parameter-file presence; it does not demonstrate binder optimization. No replacement generation job has been submitted during surface review. All caches and AF2 parameters remain in scratch. See the [Hoffman2 run record](hoffman2_pilot.md) for the job history and current outcome.

RFD4 job 1310910 completed and produced both files, but structure quality is poor. Full-target C-alpha RMSDs from 37NT are 23.264 and 21.294 A. Independent FNIIIb-domain RMSDs are 10.860 and 4.060 A, so the mismatch includes domain distortion, beyond changes in interdomain angles. Adjacent target C-alpha spacings range from 0.534–6.201 A and 0.429–5.996 A; 93 and 60 adjacent pairs, respectively, lie outside the diagnostic 3–4.5 A interval. The experimental reference ranges from 2.936–3.843 A and has three pairs just below 3 A; the interval is a diagnostic flag rather than a universal geometry rule. The severe generated spacings and domain distortions still warrant rejection. These measurements do not pinpoint the upstream cause. The experimental input sequence/coordinates and conditioning annotations were checked. Preserve the outputs as failures while auditing conditioning/decoding before further generation.

AF3 environment and official weights are installed in Hoffman2 scratch. Its assessment was canceled after the missing BindCraft arm; no AF3 prediction results exist for this cohort.

## Sources

- [MERTK 37NT](https://www.rcsb.org/structure/37NT) and [structural study](https://doi.org/10.1042/BCJ20260584).
- [GAS6 1H30](https://www.rcsb.org/structure/1H30), [AXL–GAS6 2C5D](https://www.rcsb.org/structure/2C5D), [PROS1 EGF3–4 1Z6C](https://www.rcsb.org/structure/1Z6C).
- [Predicted PROS1](https://alphafold.ebi.ac.uk/entry/P07225), [predicted GAS6](https://alphafold.ebi.ac.uk/entry/Q14393). AlphaFold DB models are licensed CC BY 4.0; cite the database and source models when reusing them.
- Reviewed UniProt sequences: human MERTK Q12866, mouse MERTK Q60805, human AXL P30530 and TYRO3 Q06418. Residue correspondence and source snapshots are included.
