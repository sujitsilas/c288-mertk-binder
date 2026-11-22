"""Can the placed Ig1 binders or an eight-residue terminal peptide reach the MERTK cut?"""
import argparse
import csv
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser, Superimposer


def chain(path, name):
    return PDBParser(QUIET=True).get_structure(path.stem, path)[0][name]


def distance(point, coords):
    return float(np.linalg.norm(coords - point, axis=-1).min())


def run(args):
    receptor = chain(args.native, "A")
    for context in args.context:
        model = chain(context, "A")
        # Place the modeled stalk in the observed membrane-proximal domain frame.
        residues = [i for i in range(386, 480) if i in receptor and i in model
                    and "CA" in receptor[i] and "CA" in model[i]]
        assert len(residues) >= 70 and model[490].resname == "PRO" and model[491].resname == "SER"
        fit = Superimposer()
        fit.set_atoms([receptor[i]["CA"] for i in residues], [model[i]["CA"] for i in residues])
        rot, shift = fit.rotran
        bond = (model[490]["C"].coord + model[491]["N"].coord) / 2
        cut = bond @ rot + shift
        for path in args.binder:
            binder = chain(path, "Z")
            endpoints = list(binder.get_residues())
            atoms = np.array([a.coord for a in binder.get_atoms() if a.element not in {"H", "D"}])
            nearest = distance(cut, atoms)
            nearest_anchor = distance(receptor[484]["CA"].coord, atoms)
            assert nearest > 0 and any(distance(receptor[i]["CA"].coord, atoms) < 12
                                       for i in range(95, 187) if i in receptor)
            # Six full CA steps from 484 to 490 plus <=3 A to the bond midpoint.
            bound = max(0, nearest_anchor - (6 * 3.9 + 3))
            yield [path.stem,context.stem,len(residues),fit.rms,nearest,
                   distance(cut,endpoints[0]["CA"].coord),distance(cut,endpoints[-1]["CA"].coord),
                   nearest_anchor,bound,3.9*(len(args.core)+1),args.core]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--context", type=Path, nargs="+", required=True)
    parser.add_argument("--binder", type=Path, nargs="+", required=True)
    parser.add_argument("--core", default="GPQGIWGQ")
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    with args.out.open("w") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["parent","stalk_model","aligned_CAs","alignment_RMSD_A","closest_binder_atom_A",
                    "N_terminal_CA_A","C_terminal_CA_A","nearest_to_observed_484_A",
                    "conservative_cut_distance_lower_bound_A","generous_core_alone_reach_A","core"])
        for row in run(args):
            w.writerow([f"{v:.3f}" if isinstance(v,float) else v for v in row])
