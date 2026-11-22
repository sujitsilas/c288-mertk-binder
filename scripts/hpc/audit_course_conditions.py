#!/usr/bin/env python3
"""Assert the installed RFD4 parser retained the intended crop and hotspot conditions."""

import argparse
import json

import rfproteina.conditions  # noqa: F401
from atomworks.io import parse
from atomworks.ml.conditions import C_SEQ
from course_condition_checks import audit_arrays
from rfproteina.conditions import C_CRDGRP, C_HOT
from rfproteina.utils.parse_presets import ANNOTATED_CIF


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input")
    parser.add_argument("output")
    parser.add_argument("--reference", required=True)
    args = parser.parse_args()
    aa = parse(args.input, config=ANNOTATED_CIF)["assemblies"]["1"][0]
    with open(args.reference) as f:
        reference = json.load(f)
    report = audit_arrays(
        chain_id=aa.chain_id,
        res_id=aa.res_id,
        res_name=aa.res_name,
        atom_name=aa.atom_name,
        coord=aa.coord,
        known=C_CRDGRP.annotation(aa, default="generate") == 1,
        sequence=C_SEQ.mask(aa),
        hot=C_HOT.mask(aa),
        reference=reference,
    )
    with open(args.output, "w") as f:
        json.dump(report, f, indent=2)
        f.write("\n")
    print(
        "COURSE_CONDITIONS_VERIFIED",
        {k: v for k, v in report.items() if k != "parsed_to_author_residue_map"},
    )


if __name__ == "__main__":
    main()
