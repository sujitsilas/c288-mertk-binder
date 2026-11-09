# Designing ligand-compatible MerTK binders: can adapted ESM3 match ProteinMPNN?

Working revision, 2026-11-09, responding to Professor Abedi's proposal feedback.
The submitted PDF remains the original proposal. This draft separates the protein
we want to build from the computational method we will test.

## What I want to build, and why

I want to develop a small, modular binder to the extracellular region of MerTK
as a building block for testing control of macrophage efferocytosis. The immediate
protein-design objective is selective MerTK binding at an accessible epitope that
permits natural-ligand engagement. A later bivalent version would test whether
receptor engagement in that geometry increases apoptotic-cell internalization
without suppressing ligand-supported clearance. Binding, receptor phosphorylation,
and improved efferocytosis are distinct outcomes; none can substitute for another.

This distinction motivates the design problem. In arthritis models, a MerTK
agonistic antibody interfered with efferocytosis despite activating the receptor.
That result does not establish what a new binder will do, but it makes preserving
ligand function a requirement rather than an optional filter.
[Waterborg et al., 2018](https://doi.org/10.3389/fimmu.2018.00742).

A compact binder is a chosen engineering platform for varying epitope, attachment
points, and valency. Its superiority to an antibody or natural-ligand stimulation
is not assumed. The eventual functional comparison must include PROS1 alone and
PROS1 plus binder, as well as an untreated condition. A candidate that binds
strongly but inhibits clearance would fail the intended biological objective.

## How I will design it

First, I will specify a MerTK extracellular construct and candidate epitope using
the available structural evidence, full ectodomain placement, glycosylation sites,
and natural-ligand access. I will compare the corresponding AXL and TYRO3 surfaces.
If no defensible patch satisfies these constraints, I will revise the protein
concept before optimizing sequences. Neither an adapter nor AF3 can correct an
inappropriate epitope or prove a productive receptor arrangement.

I will generate a common pool of candidate binder complexes using BindCraft,
then compare interface-sequence redesign on exactly the same backbones, target
constructs, and fixed-position masks:

```text
Shared target + binder complex backbones
    ├── ProteinMPNN sequence redesign ───────────┐
    └── ESM3 + LoRA + atom-adapter redesign ────┤
                                               ↓
                             Same AF3 complex-prediction protocol
                                               ↓
                     Fold/pose consistency + biological design criteria
```

The initial scope replaces the sequence-design step. AF3 is the shared structural
assessment method. The current ESM3 implementation redesigns interface residues;
it does not replace the whole structure-generation and prediction pipeline.
BindCraft's original sequences remain controls and are assessed by the same AF3
protocol. The comparison is conditional on the shared backbone pool, which already
reflects BindCraft's design and filtering choices.

## Why test an ESM adapter?

The computational question is whether a pretrained protein language model can
approach a specialized inverse-folding model for this constrained redesign task,
and whether additional target-side-chain information helps. The proposed adapter
expresses nearby target atoms in each binder residue's local frame; LoRA adapts
ESM3 while its base weights stay frozen. Training uses masked interface residues
from family-separated complexes, excluding the TAM family from model selection.

This is a hypothesis, not a requirement for making MerTK binders. Frozen ESM3 and
ESM3+LoRA distinguish adaptation from the added atomic input. ProteinMPNN is the
practical comparator. If the adapter improves over LoRA alone, I will test whether
that improvement survives a backbone-only atom control and whether it appears
with fewer training families. Merely matching ProteinMPNN would establish
feasibility; it would not establish an advantage or justify extra cost. Atomic
conditioning is also available in
[LigandMPNN](https://github.com/dauparas/LigandMPNN#22---ligand_mpnn_use_side_chain_context),
which becomes a relevant additional comparator if atomic context is claimed as
the method's advantage. Diffusion is deferred until the core question is answered.

## What would count as success?

The computational deliverable is a matched benchmark and a shortlist with explicit
pass/fail evidence for each design criterion. The primary comparison is the fraction
of sequences that retain a plausible binder fold and intended complex pose under
a common AF3 protocol. Sequence recovery is a training diagnostic. A proposed
definition of "close" is a structural pass rate no more than five percentage
points below ProteinMPNN, with paired uncertainty across independent target
families. This margin and structural cutoffs must be fixed before test evaluation;
a small, inconclusive pilot cannot establish equivalence. AF3 confidence is a
prediction-quality measure, not a binding-affinity measurement.
[AF3 output definitions](https://github.com/google-deepmind/alphafold3/blob/main/docs/output.md).

Future experimental validation would test expression and folding, MerTK binding,
AXL/TYRO3 cross-reactivity, and natural-ligand competition before testing cellular
function. The functional endpoint is increased internalization of apoptotic cells
in human macrophages, comparing monovalent and bivalent formats, a nonbinding
control, and natural-ligand stimulation. MerTK-deficient cells test receptor
dependence; phosphorylation is supporting evidence. The computational project
cannot establish improved wound healing or efferocytosis without those experiments.

Detailed criteria, the draft statistical endpoint, and implementation gaps are
recorded in [the experiment plan](../docs/experiment_plan.md).
