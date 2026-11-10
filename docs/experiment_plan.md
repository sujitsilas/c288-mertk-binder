# Design brief and comparison protocol

Reoriented on 2026-11-09 after Professor Abedi's feedback. The
[revised proposal draft](../Report/revised_project_direction.md) states the protein
concept and motivation. The working computational comparison is **ESM3+adapter
versus ProteinMPNN, with AF3 assessing both**. The shared AF3 preparation/launcher
and structural scorer are implemented and tested against output fixtures. Actual
AF3 inference has not run. See [the structural evaluation guide](af3_evaluation.md).

## Protein design comes first

The intended building block is a small MerTK-selective extracellular binder that
preserves natural-ligand access. A bivalent format is a subsequent hypothesis about
receptor engagement, not an automatic consequence of obtaining a binder. The
longer-term functional objective is greater apoptotic-cell internalization while
preserving ligand-supported efferocytosis. The course deliverable is a defensible
computational design and testing strategy; cellular activity requires experiments.

| Requirement | Why it matters | Computational evidence to collect | Decisive experimental evidence |
| --- | --- | --- | --- |
| Folded, tractable binder | A designed sequence must form a usable protein | Fold consistency, exposed hydrophobicity and aggregation flags | Expression, monodispersity and folding/stability measurements |
| MerTK binding at the selected epitope | The chosen surface must be reachable on the receptor | Intended pose and contact retention; full-ectodomain, glycan and membrane placement | Direct binding and epitope/competition measurements; cell-surface binding |
| Selectivity over AXL/TYRO3 | Family-related binding would obscure receptor-specific function | Matched off-target predictions and surface comparisons as screening evidence | Matched binding assays against all three receptors |
| Natural-ligand compatibility | A binder should not defeat the clearance mechanism it is intended to support | Overlay ligand-bound structural evidence; flag uncertain or conflicting geometries | GAS6/PROS1 competition and ligand-supported function with versus without binder |
| Productive format and valency | Strong binding and receptor activation need not improve clearance | Attachment-point, linker and ectodomain geometry assessment | Monovalent versus bivalent constructs at comparable binding-site exposure |
| MerTK-dependent efferocytosis | This is the proposed cellular benefit | No current computational metric establishes it | Apoptotic-cell internalization, distinguished from surface attachment, in control and MerTK-deficient macrophages |

Choose and document the exact construct, patch, residue numbering, structural
confidence, glycosylation assumptions, and ligand-access rationale **before**
backbone generation. The 512-token sequence-design crop is not the receptor's
physical environment: restore the candidate to a full ectodomain/membrane context
for biological screening. No patch has been selected yet. If its accessibility or
ligand compatibility cannot be justified, revisit the patch or protein concept.

The functional comparison should include ligand alone, binder alone, ligand plus
binder, a nonbinding control, and monovalent/bivalent versions. Test whether the
binder adds benefit to ligand-supported clearance, not only whether it differs
from untreated cells. Receptor phosphorylation supports the mechanism but is not
the primary functional readout. The original proposal's warning is substantive:
[Waterborg et al.](https://doi.org/10.3389/fimmu.2018.00742) observed interference
with efferocytosis from an agonistic MerTK antibody in arthritis models.

## Core computational experiment

Hold the target construct, backbone pool, interface mask, fixed residues, and
candidate budget constant. Only the sequence-design method changes. Use the same
AF3 version, weights, input preparation, seeds, number of samples, and ranking rule
for every arm. The intended comparison is conditional on a fixed backbone pool;
it is not a comparison of unrestricted end-to-end binder generators.

| Arm | Question answered | Priority |
| --- | --- | --- |
| ProteinMPNN → AF3 | How well does the established sequence-design method perform? | Core comparator |
| Frozen ESM3 → AF3 | What does the pretrained model provide without adaptation? | Core control |
| ESM3+LoRA → AF3 | Does task adaptation help without the atom adapter? | Core control |
| ESM3+LoRA+atom adapter → AF3 | Does target atomic context help and approach ProteinMPNN? | Core hypothesis |
| Original BindCraft sequence → AF3 | Does redesign help or harm the generator's starting design? | Application control |
| Adapter-only, backbone-only atoms, atom-only decoder | Is a gain attributable to LoRA, side chains, or the pretrained model? | Follow-up ablations |
| Time-conditioned masking / diffusion extensions | Does a different objective add value after the core result? | Deferred |

`configs/comparison_core.json` runs the four core sequence-design arms across three
seeds: 12 evaluation jobs and six training runs. The original eight-arm suite stays
available. See [commands and model inputs](comparisons.md). Sequence runners train
on reconstruction and export matched candidates. The separate AF3 workflow consumes
these candidates with a frozen input/sampling protocol.

ESM3 already receives backbone geometry; the adapter adds local target heavy-atom
context, not the model's first structural information. ProteinMPNN is a backbone
baseline, so a gain over it alone does not isolate an ESM-specific advantage.
[Fixed-side-chain conditioning in LigandMPNN](https://github.com/dauparas/LigandMPNN#22---ligand_mpnn_use_side_chain_context)
should be an additional comparison before claiming an advantage over atomic-context
inverse folding. It is not integrated in this repository yet.

The primary causal contrast is combined versus LoRA on matched training data and
sampling budgets. Combined versus ProteinMPNN assesses practical competitiveness.
If combined only matches ProteinMPNN while requiring more resources, report
feasibility rather than a computational advantage. If it fails to improve on LoRA,
there is no evidence to retain this adapter. A negative result still answers the
methods question, and ProteinMPNN remains usable for the protein-design goal.

## Define "close" before inspecting test results

Let q be the fraction of generated candidates meeting the same predefined
structural pass criteria, averaged within targets/families and then across
independent families. The proposed comparison is:

`delta = q(ESM3 + LoRA + adapter) - q(ProteinMPNN)`

A **draft** practical margin is five absolute percentage points: declare
noninferiority only if the lower bound of the paired 95% interval exceeds -0.05.
This is a proposed engineering tolerance, not a literature-established biological
threshold. Select and justify it before the test set is evaluated. A nonsignificant
difference or overlapping error bars does not establish equivalence. A pilot with
few independent families is exploratory and may be unable to resolve this margin.
MerTK alone is an application case, not multiple independent targets just because
many backbones or AF3 samples were generated for it.

Starting structural gates for validation calibration, **not validated hit rules**:

- Target–binder pair ipTM at least 0.8, using explicit chain identities.
- Binder Cα RMSD at most 2 Å after fitting the binder to its designed backbone.
- Binder Cα RMSD at most 3 Å after fitting the target domain to the reference target;
  this tests the binder's pose, which a binder-only fit would hide.
- No severe AF3 clash flag; inspect local interface clashes separately, since the
  global flag is not a complete interface-quality check.

Freeze atom/residue mappings, alignment domains, missing-residue handling and all
cutoffs after validation calibration and before the final test. Preserve all raw
confidence and geometry measurements. The interface epitope/contact retention,
interchain PAE, binder confidence, diversity of passing sequences, and compute per
passing candidate are secondary outcomes. These are measures of predicted
structural consistency, not measured affinity or biological function.

Select one AF3 prediction per sequence by the same predeclared AF3 ranking rule,
with equal seeds/samples for each method. Never choose by closeness to the desired
backbone after seeing outcomes. Report seed robustness separately. Count sequences
that finish prediction but fail structural gates as failures; record infrastructure
failures separately and rerun them under the same budget instead of silently
removing difficult cases. Report valid-output rate as well as structural pass rate.
AF3 describes these fields as [confidence and ranking measures](https://github.com/google-deepmind/alphafold3/blob/main/docs/output.md).

## Evaluation stages and resource budget

1. Expand and audit family-split PDB complexes. Masked native reconstruction trains
   the model and checks learning; recovery and diversity are secondary diagnostics.
   Split both partners' families, reserve TAMs, and tune no hyperparameters on test.
2. Run the four core arms. Tune each method's sampling temperature on validation
   under the same search budget; the same numeric temperature need not give the
   same diversity. Keep a common locked evaluation protocol for final reporting.
3. Run a small predeclared AF3 integration pilot, then a family-held-out structural
   benchmark. Cache the same target MSA for all candidates and use the same policy
   for binder MSAs. Use template-free prediction for the primary comparison and
   do not supply the designed complex as a template. Freeze software/database
   versions and document pretraining-overlap limitations. AF3 supports
   [explicit MSA and template settings](https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md).
4. Apply the locked protocol to shared MerTK backbone candidates. Rescore the
   original BindCraft sequence with AF3 alongside all redesigns; do not compare
   BindCraft's internal AF2 scores directly with AF3 scores. Its accepted backbone
   pool is already selected by BindCraft, so results are conditional on that pool.
5. Assess the biological criteria in the table, then rank candidates for proposed
   experimental tests. An AF3 pass never automatically satisfies those criteria.

Choose the structural cohort before seeing test confidence scores. Use all planned
candidates or a deterministic equal-sized sample per method/backbone; do not select
each method's favorable cases after prediction. Keep the original controls.
The AF3 budget is:

`targets × backbones/target × arms × design seeds × sequences/seed × AF3 seeds × samples/AF3 seed`

For example, a plumbing pilot of 3 targets × 2 backbones × 4 arms × 1 design seed
× 4 sequences × 1 AF3 seed × 5 samples requires 480 predicted structures, plus
original-sequence controls. These are not 480 independent observations. This is
an illustrative budget, not a submitted job or a power calculation. Set a feasible
budget before a full run. Measure sequence-generation and AF3 costs separately,
and report total cost including training; LoRA reduces trainable parameters but
does not remove ESM3 backbone inference cost.

Before attributing gains to useful atomic information, assess sensitivity to
plausible target-side-chain conformations: native bound side chains can provide
more favorable context than an uncertain application target. The backbone-only
control helps isolate information content but does not establish robustness.

If sample efficiency is pursued, train LoRA and combined on the same nested
family-level subsets at matched training budgets and compare both learning curves.
That follow-up can justify adaptation with limited data; parity at one data size
alone cannot establish it.

## Geometry and diffusion

[ESM3](https://github.com/Biohub/esm) already has backbone geometric reasoning and
discrete structure tokens. The new adapter supplies the additional target atomic
context proposed in the PDF. Rigid-motion invariant local coordinates preserve
orientation information that pure distances would lose. Missing or degenerate
frames and empty neighborhoods return zero conditioning. Atom inputs never
include binder side chains, including the unmasked ones.

[AlphaFold 3](https://www.nature.com/articles/s41586-024-07487-w) motivates
conditioning at both atom and residue scales and separating structural context
from a denoising task. The adapter here uses local atom pooling, not a reproduction
of AF3's Pairformer or all-atom diffusion architecture. AF3 weights are unnecessary
for training this adapter and were not downloaded.

[RFD4-Proteina](https://rfd4-proteina.org/latest/architecture/index.html) provides a
broader sequence/structure co-design reference. Its atom/latent representations
and coordinate/latent generation are substantially larger in scope than this
fixed-backbone project. The implemented extension changes the mask rate over time
and uses an absorbing-mask sequence sampler. It has no learned atomic coordinate
denoiser, no RFD4 latent autoencoder, and no flow-matching objective.

A later extension could predict continuous sequence latents conditioned on the
same geometry, or denoise binder coordinates while retaining target conditioning.
Either would require a new objective, validity checks, and independent evaluation.
Do not describe those extensions as already implemented.

## Data quality and leakage

[RCSB Search](https://search.rcsb.org/) is used for discovery;
[assembly mmCIF downloads](https://www.rcsb.org/docs/programmatic-access/file-download-services)
provide coordinates. Biological assemblies are used instead of treating arbitrary
crystal contacts as interacting pairs. Assembly annotation and proximity alone do
not establish physiological binding.

The initial pipeline uses resolved standard protein heavy atoms, highest-selected
alternate locations as handled by Bio.PDB, and full polymer sequences. It excludes
non-protein context from the learning features and skips design
positions missing N/CA/C/O (O is required by standard ProteinMPNN).
Resolution filtering applies to search results; explicit-ID downloads
bypass it. Glycans, ligands, membrane context, covalent chemistry, protonation, and
side-chain uncertainty are not modeled in this first implementation. Protein-only
polymer filtering does not mean the original complex lacks glycans or small molecules.

The split procedure joins both partners' sequence clusters, and every assembly
of the same PDB entry, before assignment. Human TAM seeds come from UniProt
[MERTK Q12866](https://www.uniprot.org/uniprotkb/Q12866/entry),
[AXL P30530](https://www.uniprot.org/uniprotkb/P30530/entry), and
[TYRO3 Q06418](https://www.uniprot.org/uniprotkb/Q06418/entry), including annotated
domains and extracellular regions to help match fragments.

30% identity with 80% bidirectional coverage is a **sequence-based approximation**
to a protein family. Remote homologs, different domain architectures, and partial
constructs may escape it. Before a scientific benchmark, audit cross-fold hits
with sensitive sequence/domain and structure searches, explicitly review TAM
ortholog/domain exclusions, and record a release-date cutoff. The seed-domain
approach can also conservatively exclude unrelated proteins with similar domains.
ESM3's pretraining overlap with public PDB structures is not eliminated by this
fine-tuning split; report it as a limitation of the frozen model baseline.

## Remaining benchmark and MerTK application work

- Expand and curate the pilot dataset; ensure enough independent validation and
  test families. Audit assemblies, missing residues, duplicate interfaces, sequence
  similarity, and domain-level leakage.
- Run the four-arm core sequence comparison first. Validate the implemented AF3
  preparation, execution and scorer on actual predictions before measuring the
  revised endpoint. Calibrate and freeze thresholds using validation data, then
  run the matched structural benchmark. Follow up with broader ablations only as
  motivated by the core result.
- Measure speed, peak CUDA memory, and numerical behavior on the lab GPUs before
  long runs. The CUDA/bf16 code path has not been tested on this laptop.
- Use the BindCraft setup/import/application runner with a selected MerTK
  extracellular construct and patch. Its source and launcher are prepared;
  generation has not been run. A full receptor sequence may exceed this pipeline's 512-token
  budget; choose and document an extracellular domain construct rather than silently
  cropping the interface. Native PDB proteins are reconstruction data, not a supply
  of de novo MerTK minibinder candidates.
- Rank candidates using separate complex prediction, AXL/TYRO3 comparisons,
  ligand-access and glycan/membrane compatibility checks. Current sequence sampling
  alone provides none of these guarantees.
- Keep bivalency, signaling, ligand competition, and macrophage internalization as
  future experimental questions, following the proposal's controls. Better native
  sequence recovery is not evidence of improved efferocytosis.
