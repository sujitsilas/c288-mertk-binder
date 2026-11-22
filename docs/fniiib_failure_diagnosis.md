# FNIIIb failure diagnosis and Ig1 priority

November 13, 2026. This is a post hoc analysis of rejected structures, not candidate nomination or evidence of measured binding. The queued Ig1 experiment remains the primary next test.

## What failed

FNIIIb job 1573763 completed 30 attempts, with a 31st interrupted by the four-hour limit. Fourteen trajectories reached structural scoring, generating 280 saved sequence proposals. Of these, 267 failed base AF2 validation; 13 reached full scoring, all from the single 130-residue trajectory `FNIIIb_F472_l130_s910883`. No variant was accepted.

| Measurement for the 13 fully scored variants | Observed range | Frozen filter |
|---|---:|---:|
| Mean interface pTM | 0.73-0.79 | >=0.50 |
| Mean normalized interface PAE | 0.27-0.33 | <=0.35 |
| Mean complex pLDDT, normalized | 0.91 | >=0.80 |
| Mean buried unsatisfied polar-atom score | 8.0-12.5 | <=4 |
| Mean shape complementarity | 0.55-0.64 | >=0.60 |
| Mean binder RMSD to designed binder, A | 0.82-1.15 | <=3.5 |

The universal final rejection is the buried-polar score. Six variants also fail shape-complementarity criteria and one fails surface hydrophobicity; these categories overlap. A plausible predicted fold and interface pose therefore coexist with an unfavorable modeled interface. The score counts buried hydrogen-bond donors/acceptors without suitable satisfaction, using the upstream ddG-style filter; it does not identify 8-12.5 individual residues or measure affinity. See the [Rosetta filter documentation](https://docs.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Filters/filter_pages/BuriedUnsatHbondsFilter).

`mpnn_fix_interface=true` fixed all 37 interface positions in the parent trajectory. Comparing the 13 full-scoring sequences confirms zero substitutions at these positions. They are correlated sequence variants of one interface, not 13 independent solutions. Additional core/surface sampling alone cannot directly replace the fixed interface chemistry. [BindCraft's upstream guidance](https://github.com/martinpacesa/BindCraft/wiki/De-novo-binder-design-with-BindCraft) describes this self-consistency step and allows interface redesign as a separate protocol choice.

## A second, independent problem: the binder occupies a glycan-bearing junction

The 13 rejected PDBs contact all three intended hotspots, V469/P471/F472, in their generated coordinates. Their receptor footprints also include numerous FNIIIa positions (293-347) as well as FNIIIb positions (382-472). This is a large interface spanning the FNIIIa/FNIIIb junction, rather than a footprint confined to the small patch used to rank the target site.

After fitting the target to the native experimental crop, all 13 structures overlap resolved carbohydrate. To test whether this was an artifact of fitting two moving domains together, each structure was also fitted separately on FNIIIa and FNIIIb:

| Alignment | Fitted target C-alpha RMSD, A | Closest binder-glycan heavy atoms, A |
|---|---:|---:|
| Whole FNIII pair | 2.65-3.59 | 0.21-0.55 |
| FNIIIa alone | 0.73-0.87 | 0.09-0.51 |
| FNIIIb alone | 0.79-0.95 | 0.31-0.66 |

With FNIIIb-only alignment, the closest glycan is **chain G, attached to N294**, for every variant. The NAG C1-Asn ND2 distance is 1.447 A, supporting the attachment. All closest distances are far below the existing 2.5 A observed-glycan screening threshold. The individual domains are reasonably preserved; the high whole-crop RMSD largely reflects their relative arrangement. The glycan conflict persists under local alignment and should not be dismissed as a global-fit artifact.

This comparison uses the resolved 37NT glycan conformation. It does not establish impossibility under every receptor/glycan conformation. It does establish that these poses do not fit the current experimental glycosylated reference. Glycans were absent during generation, so their overlap did **not** directly trigger the earlier AF2/Rosetta failures: it is an additional obstacle to native binding.

![A rejected FNIIIb design overlapping the N294 glycan](../outputs/analysis/20261113/fniiib_glycan_overlap.png)

Representative mpnn3 variant after fitting only FNIIIb. Green: rejected binder and its translucent surface. Blue: native FNIIIb. Gray: native FNIIIa. Orange spheres: native N294 glycan, extending into the binder. Magenta sticks: V469/P471/F472. The generated target is hidden; neither source coordinates nor acceptance decisions were changed. An editable [PyMOL session](../outputs/analysis/20261113/fniiib_glycan_overlap.pse) is retained.

## Why prioritize the new Ig1 patch

The screen identifies I122/F147/Y148 as three exposed preferred anchors around Q124, whereas the selected V469-centred FNIIIb patch contains F472 as its one preferred anchor. The summed target-patch SASA is approximately 1,042 versus 643 square angstroms. These observations support a stronger design hypothesis at Ig1; they are not measured energetic hotspots or a validated success predictor. The old K178/E183/V185 Ig1 pilot tested a different surface and does not answer this new question.

The Ig1 patch still requires glycan review: 31% of its screened patch members lie inside the assumed sequon exclusion radius, particularly near N114. It is also close to ligand-implicated positions, making ligand competition a meaningful test. The structural paper discusses N114 near the presumed Ig1 ligand-binding site; it does not establish that a designed binder at our patch will block either ligand. See [Starbird et al., 2026](https://pmc.ncbi.nlm.nih.gov/articles/PMC13660871/).

Recommendation: prioritize the already queued Ig1 run 1573761, retain the acceptance filters, and inspect the entire footprint of any Ig1 hit before scaling. At the November 13, 9:09 a.m. Pacific check, that job remained pending for GPU availability. Do not spend another broad tranche repeating the unchanged FNIIIb settings. A later FNIIIb rescue would first need a pose that avoids N294/glycan-bearing junctions, then independently test interface redesign; relaxing the buried-polar cutoff would not resolve the native-context problem. Thirty unsuccessful attempts do not establish that all FNIIIb epitopes are undesignable.

## Correction to the earlier glycan annotation

The reference contains modeled sugars at **10 of 14** sequons, not four. Direct attachment mapping identifies N114, N170, N207, N215, N234, N294, N316, N354, N389 and N395. This agrees with the structural paper. Four sites lack modeled sugars: N329, N336, N442 and N454. The current epitope descriptions and triage comments were corrected; frozen run inputs retain their original provenance.

## Evidence files

- `outputs/hpc/round_20261112_01/campaigns/FNIIIb_F472/bindcraft/designs/`: native statistics, sequence proposals, relaxed and rejected structures.
- `outputs/analysis/20261113/fniiib_native_context_diagnostics.json`: per-structure hashes, contacts and preliminary native-context diagnostics for 13 rejected variants and 14 scored trajectories.
- `outputs/analysis/20261113/fniiib_alignment_sensitivity.json`: domain-specific alignment checks and glycan attachment mapping.
- `outputs/analysis/20261113/fniiib_representative.json`: selected visualization source and alignment specification.

These are post hoc diagnostics and have not been used to promote rejected designs or alter the queued experiment. No additional GPU job was submitted.
