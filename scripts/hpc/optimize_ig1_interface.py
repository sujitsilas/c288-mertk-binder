"""Screen conservative affinity variants only after an independent geometry pass."""

import argparse
import json
from pathlib import Path

import numpy as np
from Bio.SeqUtils import seq3
from check_course_outputs import structure
from repair_ig1_clash import HOTSPOTS, contacts, review
from score_ig1_refolds import refold, write_json

from mertk.provenance import file_hash
from mertk.structures import read_design_chain

ALTERNATIVES = ("AILMV", "FYW", "STNQ", "DE", "KR")


def sites(path):
    model = structure(path)
    target = np.array([a.coord for a in model["A"].get_atoms() if a.element != "H"])
    seq = read_design_chain(path, "B")[0].sequence
    found = []
    for (number, aa), residue in zip(enumerate(seq, 1), model["B"], strict=True):
        if number in {15, 32, 83}:
            continue
        xyz = np.array(
            [a.coord for a in residue if a.element != "H" and a.name not in {"N", "CA", "C", "O"}]
        )
        choices = next((group.replace(aa, "") for group in ALTERNATIVES if aa in group), "")
        if not len(xyz) or not choices:
            continue
        distance = float(np.linalg.norm(xyz[:, None] - target[None, :], axis=-1).min())
        if 2.4 <= distance <= 4.5:
            found.append((distance, number, aa, choices[:2]))
    return sorted(found)[:3]


def mutate_repack(pose, position, identity, pr, score):
    changed = pose.clone()
    index = changed.pdb_info().pdb2pose("B", position)
    if identity:
        pr.rosetta.protocols.simple_moves.MutateResidue(index, identity).apply(changed)
    movemap = pr.rosetta.core.kinematics.MoveMap()
    movemap.set_chi(False)
    movemap.set_bb(False)
    movemap.set_jump(False)
    for r in range(
        max(changed.chain_begin(2), index - 3), min(changed.total_residue(), index + 3) + 1
    ):
        movemap.set_chi(r, True)
    relax = pr.rosetta.protocols.relax.FastRelax(score, 1)
    relax.set_movemap(movemap)
    relax.set_movemap_disables_packing_of_fixed_chi_positions(True)
    relax.max_iter(100)
    relax.apply(changed)
    return changed


def score_interface(pose, pr, score):
    mover = pr.rosetta.protocols.analysis.InterfaceAnalyzerMover()
    partners = getattr(pr.rosetta.core.pose, "DockingPartners", None)
    mover.set_interface(partners.docking_partners_from_string("A_B") if partners else "A_B")
    mover.set_scorefunction(score)
    mover.set_compute_interface_energy(True)
    mover.set_pack_separated(True)
    mover.apply(pose)
    return float(mover.get_interface_dG())


def main(suite, project):
    screens = json.loads((suite / "independent.json").read_text())
    result = suite / "affinity_screen.json"
    eligible = {name: row for name, row in screens.items() if row["pose_pass"]}
    if not eligible:
        write_json(result, {"status": "Skipped. No independently refolded pose passed geometry."})
        return
    import pyrosetta as pr

    pr.init("-mute all -constant_seed -jran 112026")
    key = min(eligible, key=lambda name: eligible[name]["binder_pose_ca_rmsd_A"])
    refolded = suite / "refolds" / key / "refold.pdb"
    if file_hash(refolded) != eligible[key]["sha256"]:
        raise ValueError("Independent pose differs from the scored file")
    directory = suite / "affinity"
    directory.mkdir(exist_ok=True)
    # Repack against the experimental receptor in the matched alignment.
    from relax_ig1_interface import placed_complex

    start = directory / "reference_pose.pdb"
    reference = suite / "target.pdb"
    placed_complex(refolded, reference, start)
    baseline = {n for n in contacts(start, reference) if 93 <= n <= 185}
    pose = pr.pose_from_pdb(str(start))
    score = pr.get_fa_scorefxn()
    scores = []
    for _, number, wt, alternatives in sites(start):
        control = score_interface(mutate_repack(pose, number, None, pr, score), pr, score)
        for aa in alternatives:
            name = f"{wt}{number}{aa}"
            path = directory / f"{name}.pdb"
            mutated = mutate_repack(pose, number, seq3(aa).upper(), pr, score)
            mutated.dump_pdb(str(path))
            row = review(path, suite, start, baseline, max_rmsd=1.0)
            ddg = score_interface(mutated, pr, score) - control
            scores.append(dict(id=name, relative_interface_energy_REU=ddg, geometry=row))
            print(name, ddg, row["pose_pass"], flush=True)
    report = dict(status="Scored conservative single substitutions", parent=key, variants=scores)
    write_json(result, report)
    improved = [
        row
        for row in scores
        if row["geometry"]["pose_pass"] and row["relative_interface_energy_REU"] < 0
    ]
    if improved:
        top = min(improved, key=lambda row: row["relative_interface_energy_REU"])
        path = directory / f"{top['id']}.pdb"
        verify = directory / f"{top['id']}_refold"
        verify.mkdir(exist_ok=True)
        refold(
            project,
            reference,
            suite / "37NT.cif",
            read_design_chain(reference, "A")[0],
            read_design_chain(path, "B")[0],
            verify,
            HOTSPOTS,
        )
        report["top_variant"] = top["id"]
        report["top_refold"] = review(verify / "refold.pdb", suite, start, baseline, max_rmsd=2.5)
        write_json(result, report)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("suite", type=Path)
    parser.add_argument("project", type=Path)
    args = parser.parse_args()
    main(args.suite, args.project)
