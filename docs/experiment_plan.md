# Proposal review and implementation decisions

The proposal's strongest first result is a controlled benchmark of target
side-chain conditioning for interface reconstruction. It separates a measurable
computational question from the more uncertain MerTK agonism/efferocytosis
application. Keep these as separate claims in the report.

## First experiment

| Arm | Frozen ESM3 | Train LoRA | Train target-atom adapter | Time-conditioned masking |
| --- | --- | --- | --- | --- |
| ESM3 baseline | yes | no | no | no |
| LoRA | yes | yes | no | no |
| Geometry | yes | no | yes | no |
| Combined | yes | yes | yes | no |
| Combined, backbone-only target atoms | yes | yes | yes, restricted atoms | no |
| Denoising extension | yes | yes | yes | yes |
| Atom-only decoder | absent | no | yes + small decoder | no |
| ProteinMPNN | absent | no | no; pretrained backbone model | no |

All ESM3 arms receive the same target sequence, binder/target backbone coordinates,
backbone structure tokens, and design masks. ProteinMPNN receives N/CA/C/O backbone
coordinates and fixed-position sequences; the small decoder receives only the
local target-atom adapter features. All arms use the same complexes and eligible
design positions. See [the comparison guide](comparisons.md) for the full input
table, three-seed runner, and separate BindCraft application comparison.

Training-time validation measures single-pass full-interface recovery;
cross-entropy uses the full ESM3 vocabulary, while recovery/generation are limited
to the 20 standard amino acids. Training masks a subset of the interface, so
reconstruction difficulty differs from the full-mask evaluation by design.
The cross-model comparison samples from a fully masked interface and reports
sample-averaged recovery and diversity. Its sampling protocol differs from
single-pass validation and must be labeled separately in results.

Run at least three seeds and compare matched family subsets for learning curves.
Aggregate confidence intervals over independent family components, not individual
residues or duplicated assembly copies. The current trainer samples examples
uniformly, so repeated crystal copies can overweight an entry; deduplicate or
add component-balanced sampling before scaling the benchmark. Never tune on the
test or TAM application splits.

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
- Run the implemented ProteinMPNN, atom-input small-decoder, backbone-only
  target-atom, and ESM3 ablation suite on the expanded dataset. Plot learning curves
  from each run's metrics and add independent complex structure validation.
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
