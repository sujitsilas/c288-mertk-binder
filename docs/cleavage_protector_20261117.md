# Dedicated MERTK cleavage-protector pilot — November 17, 2026

The previous `protector_triage` job only screened the shared Ig library and found no eligible parent. It did not generate a cleavage protector. Following the user's request to triage and run this separate aim, the new campaign designs monovalent binders against the **MERTK cleavage-site stalk**, independently of Ig1/Ig2 generation. Direct inhibition of ADAM17 is not this campaign's mechanism.

The reported human cleavage bond is **P490–S491**. MERTK-binding nanobodies have inhibited shedding without directly inhibiting ADAM17 or blocking GAS6 binding, but their steric versus allosteric mechanism remains unresolved. This supports the biological hypothesis, not our particular epitope or sequence. [Primary study](https://pmc.ncbi.nlm.nih.gov/articles/PMC13378718/).

## Target triage

- **Hotspots:** Y488, P490 and S491, within stalk A485–A498. The previously failed FNIIIb F472 patch remains retired.
- **Design context:** A386–A526, retaining the adjacent domain and transmembrane helix as context. Full A93–A526 and experimental 37NT protein context are used in structural review.
- **Source:** [AFDB Q12866 v6](https://alphafold.ebi.ac.uk/entry/Q12866), sequence-checked against the frozen canonical human reference. Experimental 37NT ends at G484 and does not observe the cleavage bond.
- **Uncertainty:** stalk pLDDT is about 33–40. Mean predicted alignment error from stalk to the adjacent domain is about 22.9 Å and to the transmembrane helix about 24.3 Å. The pose is a hypothesis, not a solved target structure.
- **Conformations:** original model and a +25° perturbation of P490 psi receive separate generation runs. A −25° conformer is reserved for sensitivity screening. Bond lengths and backbone continuity are preserved. These are geometric sensitivity models, not a physical equilibrium ensemble.
- **Accessibility probe:** 538/1024, 550/1024 and 522/1024 directions on a 5 Å sphere around the cleavage bond pass the declared protein/glycan/sequon/membrane exclusions in the three respective models. This tests contact-point access, not accommodation of a whole binder or ADAM17 exclusion.
- **During generation:** locally transported resolved glycans at N389/N395, conservative sequon exclusion at N389/N395/N442/N454, and a modeled membrane plane with a 5 Å extracellular margin. No lipids or glycans are cofolded. Full mammalian glycan accessibility remains unproven.

The earlier Ig glycan-loss runtime passed; this changed target and membrane term require their own CPU and full GPU forward/backward preflights. CPU tests passed for both target inputs, including finite, repulsive membrane gradients. BindCraft's default final acceptance filters and pinned upstream source are unchanged.

## Submitted jobs

| Stage | Job | Resource / cap |
|---|---:|---|
| Input and CPU gradient validation | 1896122 | CPU, 10 min; completed |
| Both-target GPU integration preflight | 1896133 | 1 GPU, 10 min |
| Original stalk conformation generation | 1896138 | 1 GPU, 55 min |
| Original-conformation QC | 1896152 | CPU, 15 min |
| +25° stalk conformation generation | 1896184 | 1 GPU, 55 min |
| +25° conformation QC | 1896188 | CPU, 15 min |
| Protector parent selection | 1896201 | CPU, 10 min |
| Conditional parent AF3 | 1896274 | 1 GPU, 20 min; held until selection |
| Final report | 1896276 | CPU, 15 min |

Generation permits **A100/H200/V100**; AF3 permits **A100 80 GB/H200/RTX PRO 6000 Blackwell**. The account currently routes to campus. Requests use four CPU cores; no GPU is allocated to CPU-only stages. AF3 is released only if there is an eligible parent; otherwise the selection stage cancels the held job before allocation.

Live accounting at submission: **21.4433 GPU-hours** already consumed, plus **2.3333 hours maximum** for this suite, totaling **23.7767/24 hours**. No increased budget is assumed. Generation is bounded by scheduler walltime; interruptions and actual completed attempts will be reported.

## What constitutes a result

The generator searches 55–90 residue binders. Accepted poses must contact the cleavage region, place at least half of contacted target residues within A485–A498, avoid transmembrane contacts, and pass the three conformer/context checks. Model and restored-protein clashes, local observed-glycan overlap, nearby sequons and membrane approach are screened separately. These are provisional model gates; no candidate is automatically qualified as a protective molecule.

At most one passing parent proceeds to three-seed AF3 (101/202/303, one sample each, ten recycles, no templates, single-sequence MSAs), using MERTK A93–A526. Previous TULP2/mCherry predictions used A93–A484 and are displayed only as **historical, unmatched references**; their scores do not calibrate this changed target construct.

Experimental advancement requires native MERTK binding, reduced soluble MERTK alongside retained intact surface receptor, preserved GAS6/PROS1 function and MERTK-dependent apoptotic-cell internalization. Controls must distinguish assay interference, receptor depletion and direct ADAM17 inhibition. Binding confidence and modeled cleavage-site coverage do not establish protection.

## Records and review

- [Frozen suite and job commands](../outputs/hpc/round_20261117_protector_01/manifest.json)
- [Input geometry and confidence](../outputs/hpc/round_20261117_protector_01/references/triage.json)
- [Accessibility probe](../outputs/design/cleavage_protector_20261117/accessibility_probe.json)
- [Protector and Ig pilot gallery rebuild instructions](prediction_review.md#rebuild-the-pilot-galleries)

The watcher synchronizes logs/results every two minutes while the Mac is awake and connected. Slurm runs independently. The protector viewer labels its receptor and glycan placement as modeled and only displays binder poses after actual output files exist.
