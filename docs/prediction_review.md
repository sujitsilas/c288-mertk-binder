# MERTK prediction viewer and independent comparators

The [final project review](../outputs/review/mertk_protector_predictions.html) shows the generation and repair results for both aims alongside the GAS6 and mCherry Boltz-2 comparisons. The earlier pilot galleries used TULP2 and mCherry as AF3 comparators and follow the side-by-side structure-review approach in Gokul's presentation.

## Earlier pilot gallery methods

The left panel displays every Ig1/Ig2 Accepted PDB with a matching BindCraft final-statistics row. It preserves original geometry in the predicted-complex view and offers a separately labeled overlay onto the experimental MERTK ectodomain, A93–A484. The latter retains resolved glycans and adjacent domains. It does not reconstruct the missing stalk, membrane or intracellular domain. Choose Ig1, Ig2 or both Ig domains for alignment; reported distances and RMSDs follow the selected fit. Use the original view to check whether an overlay is hiding target deformation.

The page supports rotation, zoom, interface highlighting, glycan visibility, a predicted-target overlay, PNG export and partner-coordinate export. Display coordinates omit hydrogens for clarity; original files remain unchanged. Native numbering is recovered by an exact, unique match to human MERTK sequence, not by trusting generator numbering. New filter passes remain provisional; the viewer cannot establish binding, efferocytosis or shedding protection.

| Comparator | Why included | Interpretation limit |
|---|---|---|
| Human TULP2, UniProt O00295, 520 aa | Tulp2 was negative in published Mer-Fc binding experiments. | That experimental context does not establish the exact human sequence/construct used here as a nonbinder. Low confidence in disordered regions is not proof of nonbinding. |
| mCherry precursor, 236 aa | A folded protein unrelated to the MERTK ligand-design objective; experimental structure 2H5Q documents its fold. | It is a putative unrelated-protein comparator, not an experimentally established MERTK nonbinder. The standard-amino-acid AF3 input does not model chromophore maturation. |

The Tulp2 evidence comes from [Caberoy et al., 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC3020645/). mCherry structure and sequence come from [PDB 2H5Q](https://www.rcsb.org/structure/2H5Q). The PDB's CH6 chromophore represents precursor residues MET–TYR–GLY, so the input expands that component to MYG using the deposited chemical-component definition; it does not replace it with one X or delete it.

The [control bundle](../outputs/controls/mertk_comparators_20261113/manifest.json) contains two AF3 input files. It uses the complete observed MERTK sequence, three seeds (101/202/303), one sample per seed, no templates and explicit single-sequence MSAs. The November 14 campaign returned three predictions per control. [AF3 results manifest](../outputs/hpc/round_20261114_ig_guard_01/comparisons/results_manifest.json)

After a candidate passes native-context review and is nominated, predict it alongside the controls under the same target construct, engine, MSA policy, seeds and sample budget. The current BindCraft AF2 means cannot be ranked directly against future control AF3 scores. Inspect all control seeds; a high-scoring negative is evidence about model specificity, not automatically a new ligand. Recheck the shared GPU budget before release; the existing AF3 reserve also covers controls. These inputs are not compatible with the legacy candidate-cohort launcher unchanged, because controls lack designed-reference poses and are not binder nominees.

For an actual control result, add one entry per seed to that case's `predictions` list in a separately recorded results manifest: `pose` (PDB or mmCIF path relative to the manifest), `input_sha256`, `seed`, and optional `iptm`/`partner_plddt`. Point a review config at that results manifest. The builder verifies input hashes and exact target/partner sequences. Preserve the prepared manifest and real raw prediction files; do not insert a reference structure or hypothetical docked pose as a complex result.

## Rebuild the pilot galleries

The previous gallery exports are archived. Run these commands from the repository root to rebuild them from the retained pilot inputs and predictions.

```bash
.venv/bin/python scripts/build_prediction_review.py outputs/hpc/round_20261114_ig_guard_01
.venv/bin/python scripts/build_protector_review.py outputs/hpc/round_20261117_protector_01
```

The builders write the pilot HTML and provenance JSON under `outputs/review/`. The monitor refreshes the Ig gallery when it synchronizes an active campaign.
