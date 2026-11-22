#!/usr/bin/env python3
"""Build a model-labeled stalk review; historical controls are explicitly unmatched."""

import argparse
import copy
import json
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFIO, MMCIFParser, PDBParser

from mertk.glycan_guard import SUGARS, frames
from mertk.prediction_review import build_review


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("suite", type=Path)
    args = parser.parse_args()
    suite = args.suite.resolve()
    reference = PDBParser(QUIET=True).get_structure(
        "model", suite / "references/base_full_context.pdb"
    )[0]
    native = MMCIFParser(QUIET=True).get_structure("native", suite / "references/37NT.cif")[0]
    for chain in native:
        sugars = [r for r in chain if r.resname in SUGARS]
        if not sugars:
            continue
        candidates = sorted(
            (float(np.linalg.norm(s["C1"].coord - r["ND2"].coord)), r.id[1])
            for s in sugars
            if s.resname == "NAG" and "C1" in s
            for r in native["A"]
            if r.resname == "ASN" and "ND2" in r
        )
        if not candidates or candidates[0][0] > 2.2:
            raise ValueError("Unknown glycan attachment")
        n = candidates[0][1]
        original = np.array([native["A"][n][a].coord for a in ["N", "CA", "C"]])
        target = np.array([reference["A"][n][a].coord for a in ["N", "CA", "C"]])
        copied = copy.deepcopy(chain)
        for residue in list(copied):
            if residue.resname not in SUGARS:
                copied.detach_child(residue.id)
        copied.id = f"G{len(list(reference))}"
        for atom in copied.get_atoms():
            atom.coord = (atom.coord - original[1]) @ frames(original).T @ frames(target) + target[
                1
            ]
        reference.add(copied)
    out = suite / "review"
    out.mkdir(exist_ok=True)
    writer = MMCIFIO()
    writer.set_structure(reference)
    writer.save(str(out / "modeled_context.cif"))
    config = json.loads(Path("configs/prediction_review.json").read_text())
    config.update(
        reference=str(out / "modeled_context.cif"),
        canonical_sequence=str(suite / "references/canonical.json"),
        output="outputs/review/mertk_protector_pilot.html",
        additional_suites=[],
        queue_suites=[],
        rfd4_suites=[],
        parent_refolds={},
        include_protector=True,
        modeled_reference=True,
        reference_note="Predicted MERTK A93–A526 from AFDB v6 with observed glycan shapes transported into local Asn frames. Stalk pLDDT ~33–40. Stalk, membrane orientation and glycan placement are modeled, not experimental.",
        control_note="Historical TULP2/mCherry predictions used MERTK A93–A484. They are reference comparisons, not matched controls for the stalk/TM construct. AF3 scores cannot establish inhibition of cleavage.",
    )
    path = out / "viewer_config.json"
    path.write_text(json.dumps(config, indent=2) + "\n")
    print(json.dumps(build_review(suite, path)))


if __name__ == "__main__":
    main()
