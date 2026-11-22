# Comparing established design tools on MERTK

Tool and source review: November 10, 2026. The biological question and experimental criteria are in the [proposal](repository_archive.md#proposal-drafts) and [validation plan](experiment_plan.md).

## Tool roles and practical choice

| Tool | Role | Decision for this project |
|---|---|---|
| BindCraft | Generates binders using AF2-guided design and its native sequence-design/filtering pipeline | Initial route, supported by published experimental binder results; preserve its original sequences |
| RFD4-Proteina | Generates protein sequence and structure with target conditioning | Independent comparator using the same experimental MERTK construct; selectivity and MERTK activity remain unproven |
| LigandMPNN | Sequence design on an existing structure, optionally using fixed side-chain context | Optional paired refinement of selected parent designs; keep target sequence fixed |
| AF3 | Predicts a supplied sequence complex and reports structural confidence | Common assessment after generation; not a substitute for a binder generator or for experimental receptor coordinates |

Sources: [BindCraft](https://doi.org/10.1038/s41586-025-09429-6), [RFD4-Proteina](https://github.com/RosettaCommons/RFD4-Proteina), [LigandMPNN](https://github.com/dauparas/LigandMPNN), [AF3 inputs](https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md).

BindCraft is the practical starting choice, not a claim of biological superiority. There is no measured MERTK result here that establishes a winner. AF3 plus LigandMPNN alone does not replace the scaffold-generation problem. RFD4-Proteina already supplies sequences, so mandatory redesign would obscure what its native pipeline produces.

## Shared conditions and honest denominators

Use the same experimental receptor file, sequence, reviewed epitope class, binder-length range (initially 65–120 residues), and final scoring protocol. Give each generator a prespecified GPU-hour allowance on the same GPU type and account for search/refolding costs. A proposed pilot allowance is 8 GPU-hours per method per epitope, across three recorded design seeds; revise this after a runtime-only pilot and before comparative generation. Report actual cost and failed attempts if either method cannot complete the allowance.

Archive total trajectories, generated structures, native-filter rejects, numerical failures, duplicates, completed AF3 assessments, nominations, expressed proteins and experimental hits. Report native yield and shared-filter yield separately. Retain original outputs; never compare only each tool's accepted designs and call that an overall success rate. Freeze clustering/diversity and nomination rules before examining experimental outcomes.

AF3 confidence, RMSD and model-native metrics have different meanings. Do not compare their raw values across generators as if they were affinities. The structural scorer reports only the imported cohort's outcome; total-generation yield and GPU-hour efficiency require the complete upstream attempt ledger. The code permits unequal cohort sizes but requires every arm to cover the same campaign/epitope strata. Any departures from the planned budget belong in the report.

For LigandMPNN, run a separate paired follow-up on identical parent backbones from both generators, with identical binder design positions, seeds and sequence counts. Retain `bindcraft`, `rfd4_proteina`, `bindcraft_ligandmpnn` and `rfd4_proteina_ligandmpnn` arm labels. Original versus refined sequences are paired by parent; unrelated de novo backbones are not paired observations. Do not treat multiple seeds on MERTK as independent receptor families.

## Source setup and generation

Run repository utilities from the project root in the local environment. Upstream engines need separate supported environments; the Mac environment prepares inputs but does not run the CUDA pipelines.

```bash
python scripts/setup_external.py --tool bindcraft
python scripts/setup_external.py --tool rfd4_proteina
python scripts/setup_external.py --tool ligandmpnn
python scripts/setup_external.py --tool alphafold3
```

These commands clone pinned source only. Follow each pinned checkout's installation instructions for its runtime and weights; record resolved package versions and weight checksums. BindCraft has a project helper, `bash scripts/install_bindcraft.sh`, on Linux with CUDA and Conda.

For BindCraft, use `scripts/prepare_bindcraft.py` with the reviewed target PDB and author-numbered hotspots. It validates chain/residue availability and preserves the target map, settings, filters and source revision. Its requested final-design count is not a compute budget: use the externally allocated runtime and retain all trajectory logs when it expires. Do not choose hotspot numbers simply because the example syntax accepts them.

For RFD4-Proteina, use its native condition-spec interface. Create a reviewed campaign YAML with `input` set to the experimental mmCIF and a target contig corresponding to the audited construct. The published backbone-coordinate and sequence conditions are `C_CRDGRP` and `C_SEQ`. Use the [native input documentation](https://rfd4-proteina.org/latest/how-to/input_prep/yaml-json.html) and [specificity example](https://github.com/RosettaCommons/RFD4-Proteina/blob/main/docs/examples/ppi_specificity_tutorial.ipynb) to encode the intended patch; do not silently substitute a predicted receptor.

After installing the pinned RFD4 environment and preparing `campaign.yaml`, its native CLI is:

```bash
rfd4-generate inputs=campaign.yaml num_replicates=2 \
  out_dir=./runs/mertk_rfd4 '~metrics'
```

This is a small generation-only pilot with upstream metrics disabled, not the complete budgeted comparison. Record the actual seed configuration and resolved weight revision. The campaign must choose its search/refolding policy explicitly before a comparative run. Common AF3 scoring remains separate. Do not install unrelated optional scoring plugins.

For LigandMPNN, install its runtime and model parameters separately. The following upstream command illustrates a binder-only refinement, **assuming the inspected input uses target A and binder B**:

```bash
python external/LigandMPNN/run.py \
  --model_type ligand_mpnn \
  --checkpoint_ligand_mpnn external/LigandMPNN/model_params/ligandmpnn_v_32_010_25.pt \
  --pdb_path outputs/selected_parent.pdb \
  --out_folder outputs/ligandmpnn_parent \
  --chains_to_design B --ligand_mpnn_use_side_chain_context 1 \
  --seed 42 --batch_size 1 --number_of_batches 4
```

Run this with the LigandMPNN environment's Python. To restrict refinement to selected binder positions, use its documented `--redesigned_residues` and preserve those masks for each parent comparison. Fixed side-chain context can inform sequence choice; it does not guarantee correct glycans or a productive MERTK interface. Inspect the emitted sequences, chain mappings and structures before import. Do not rename an unchanged source PDB as a refined design without inserting the actual new sequence.

## Import, provenance and common assessment

Write a per-run JSON record containing `tool_revision` (required), weight identifiers/checksums, environment, resolved settings, seeds, target/epitope mapping, attempts/rejections/failures, GPU model and elapsed GPU time. The importer carries this record unchanged; it does not verify claims about external execution. Do not use a made-up seed for an unrecorded run; recover the generator log first.

```bash
python scripts/import_bindcraft.py \
  --input outputs/bindcraft_pilot/designs/Accepted \
  --output data/candidates/bindcraft \
  --campaign-id mertk_pilot --epitope-id reviewed_patch_1 \
  --target-reference data/references/37NT/37NT.cif \
  --run-record outputs/bindcraft_run.json --seed 42

python scripts/import_candidates.py \
  --input outputs/rfd4_selected_complexes --output data/candidates/rfd4 \
  --arm rfd4_proteina --campaign-id mertk_pilot --epitope-id reviewed_patch_1 \
  --target-reference data/references/37NT/37NT.cif \
  --run-record outputs/rfd4_run.json --seed 42 \
  --target-chain A --binder-chain B

python scripts/merge_candidates.py \
  --inputs data/candidates/bindcraft/candidates.jsonl data/candidates/rfd4/candidates.jsonl \
  --output data/candidates/pilot.jsonl
```

These examples require actual completed runs and inspected chain IDs. RFD4 may assign different chains; explicitly set them. Import each seed/epitope batch separately. The importer accepts PDB or `.cif`, requires complete observed protein backbones, rejects gaps instead of inventing coordinates, and stores sequence/mapping/structure hashes. Missing residues require an explicit, documented construct decision first.

Proceed to the [AF3 guide](af3_evaluation.md). Glycan/ligand/membrane review and AXL/TYRO3 negative-target assessment are additional required steps; the current two-chain scorer does not perform them. The Annexin A5–PROS1 arm is a biological comparator with a different topology, not another miniprotein generator arm to pool into these structural statistics.
