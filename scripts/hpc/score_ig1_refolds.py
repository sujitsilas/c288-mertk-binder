"""Refold a fixed Ig1 binder sequence using only the native receptor as an AF2 template."""

import argparse
import json
from pathlib import Path

import numpy as np
from check_course_outputs import analyze

from mertk.provenance import file_hash
from mertk.structural import rigid_fit, rmsd
from mertk.structures import read_design_chain

MODEL = "model_1_multimer_v3"
SEED = 11192026


def write_json(path, value):
    pending = path.with_suffix(".tmp")
    pending.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    pending.replace(path)


def rosetta_energy(path):
    import pyrosetta as pr

    pr.init("-mute all -constant_seed -jran 111926")
    analyzer = pr.rosetta.protocols.analysis.InterfaceAnalyzerMover()
    partners = getattr(pr.rosetta.core.pose, "DockingPartners", None)
    analyzer.set_interface(partners.docking_partners_from_string("A_B") if partners else "A_B")
    analyzer.set_scorefunction(pr.get_fa_scorefxn())
    analyzer.set_compute_interface_energy(True)
    analyzer.set_pack_separated(True)
    analyzer.apply(pr.pose_from_pdb(str(path)))
    return float(analyzer.get_interface_dG())


def parent_complex(target, ligand, output):
    blocks = []
    for source, chain in [(target, "A"), (ligand, "B")]:
        blocks.extend(
            line[:21] + chain + line[22:]
            for line in source.read_text().splitlines()
            if line.startswith("ATOM  ")
        )
        blocks.append("TER")
    output.write_text("\n".join([*blocks, "END", ""]))


def refold(project, target, context, design_target, design_binder, directory, hotspots):
    from colabdesign import mk_afdesign_model

    sequence = design_binder.sequence
    model = mk_afdesign_model(
        protocol="binder",
        use_multimer=True,
        model_names=[MODEL],
        data_dir=str(project / "external/BindCraft"),
        use_bfloat16=False,
    )
    model.prep_inputs(pdb_filename=str(target), target_chain="A", binder_len=len(sequence))
    print(f"Starting {MODEL} on {len(sequence)} binder residues", flush=True)
    model.predict(seq=sequence, seed=SEED, models=[0], num_recycles=2, verbose=False)
    predicted = directory / "refold.pdb"
    model.save_pdb(str(predicted), get_best=False)
    geom = analyze(predicted, target, context_path=context, hotspots=hotspots)
    if (geom["target_chain"], geom["binder_chain"], geom["binder_sequence"]) != (
        "A",
        "B",
        sequence,
    ):
        raise ValueError("AF2 output chain identity or binder sequence mismatch")
    receptor = read_design_chain(predicted, "A")[0].backbone[:, 1]
    binder = read_design_chain(predicted, "B")[0].backbone[:, 1]
    rotation, translation = rigid_fit(receptor, design_target.backbone[:, 1])
    n = len(receptor)
    pae = np.asarray(model.aux["pae"])
    if pae.shape != (n + len(sequence), n + len(sequence)):
        raise ValueError(f"Unexpected PAE dimensions {pae.shape}")
    interchain = (pae[:n, n:] + pae[n:, :n].T) / 2
    contacts = np.linalg.norm(receptor[:, None] - binder[None, :], axis=-1) <= 10
    metrics = dict(
        model=MODEL,
        seed=SEED,
        recycles=2,
        receptor_template_sha256=file_hash(target),
        ipTM=float(model.aux["i_ptm"]),
        mean_interchain_pae_A=float(interchain.mean()),
        contact_pairs_10A=int(contacts.sum()),
        mean_contact_pae_A=float(interchain[contacts].mean()) if contacts.any() else None,
        target_aligned_binder_ca_rmsd_A=rmsd(
            binder @ rotation + translation, design_binder.backbone[:, 1]
        ),
        native_context=geom,
    )
    write_json(directory / "af2_metrics.json", metrics)
    metrics["rosetta_interface_dG_REU"] = rosetta_energy(predicted)
    write_json(directory / "metrics.json", metrics)
    print(f"Completed {directory.name}: ipTM={metrics['ipTM']:.3f}", flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("phase", choices=["parents", "repairs"])
    p.add_argument("key", choices=["Ig1_mpnn12", "Ig1_mpnn8"])
    for name in ("project", "suite", "repair"):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    target = a.repair / "inputs/target.pdb"
    context = a.repair / "inputs/37NT.cif"
    manifest = json.loads((a.repair / "manifest.json").read_text())
    hotspots = manifest["experiments"][a.key]["hotspots"]
    out = a.suite / "results" / a.phase / a.key
    out.mkdir(parents=True, exist_ok=True)
    if a.phase == "parents":
        ligand = a.repair / "inputs" / f"{a.key}_binder.pdb"
        for path in (target, ligand, context):
            if file_hash(path) != manifest["input_sha256"][str(path.relative_to(a.repair))]:
                raise ValueError(f"Frozen input mismatch: {path}")
        native = out / "placed_parent.pdb"
        parent_complex(target, ligand, native)
        baseline = analyze(native, target, context_path=context, hotspots=hotspots)
        baseline["rosetta_interface_dG_REU"] = rosetta_energy(native)
        write_json(out / "placed_parent.json", baseline)
        print("Placed parent metrics saved", flush=True)
        refold(
            a.project,
            target,
            context,
            read_design_chain(target, "A")[0],
            read_design_chain(ligand, "Z")[0],
            out,
            hotspots,
        )
        return
    records = json.loads((a.repair / "quality.json").read_text())["records"][a.key]
    passing = [r for r in records if r["preliminary_geometry_pass"]]
    write_json(
        out / "selection.json",
        dict(
            generated=len(records),
            eligible=len(passing),
            selected_sha256=[r["sha256"] for r in passing],
        ),
    )
    for row in passing:
        path = Path(row["path"])
        if file_hash(path) != row["sha256"]:
            raise ValueError(f"Repair changed since geometry screen: {path}")
        candidate = out / row["sha256"]
        candidate.mkdir(exist_ok=True)
        refold(
            a.project,
            target,
            context,
            read_design_chain(path, row["target_chain"])[0],
            read_design_chain(path, row["binder_chain"])[0],
            candidate,
            hotspots,
        )
    print(f"{len(passing)} eligible of {len(records)} generated", flush=True)


if __name__ == "__main__":
    main()
