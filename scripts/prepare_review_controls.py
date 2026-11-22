#!/usr/bin/env python3
"""Prepare different, source-checked negative comparators; never fabricate control poses."""

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

from mertk.provenance import file_hash
from mertk.structures import AA, read_design_chain


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("Preserve prior control inputs; choose a fresh directory")
    source = Path("data/references/control_sequences")
    tulp2 = json.loads((source / "O00295.json").read_text())
    assert tulp2["primaryAccession"] == "O00295"
    entity = json.loads((source / "2H5Q_entity1.json").read_text())
    component = json.loads((source / "CH6.json").read_text())["chem_comp"]
    assert component["mon_nstd_parent_comp_id"] == ["MET", "TYR", "GLY"]
    assert component["one_letter_code"] == "MYG"
    # The mature chromophore encodes three precursor residues, not one X residue.
    precursor = entity["entity_poly"]["pdbx_seq_one_letter_code"].replace("(CH6)", "MYG")
    assert len(precursor) == 236 and set(precursor) <= set(AA)
    target, _ = read_design_chain("data/references/37NT/target.pdb", "A")
    cases = [
        {
            "id": "mertk_tulp2",
            "name": "TULP2",
            "role": "Reported negative comparator",
            "sequence": tulp2["sequence"]["value"],
            "accession": "O00295",
            "evidence_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC3020645/",
            "sequence_url": "https://www.uniprot.org/uniprotkb/O00295/entry",
            "caveat": "Tulp2 was negative in reported Mer-Fc binding assays. This human sequence and our receptor construct require their own validation; disordered regions may lower prediction confidence.",
        },
        {
            "id": "mertk_mcherry",
            "name": "mCherry",
            "role": "Unrelated folded protein",
            "sequence": precursor,
            "accession": "PDB 2H5Q, precursor sequence",
            "evidence_url": "https://www.rcsb.org/structure/2H5Q",
            "sequence_url": "https://data.rcsb.org/rest/v1/core/polymer_entity/2H5Q/1",
            "caveat": "Unrelated-protein comparator, not a proven MERTK nonbinder. AF3 input uses the standard-amino-acid precursor (CH6 expanded to MYG); chromophore maturation is not modeled.",
        },
    ]
    (args.output / "inputs").mkdir(parents=True)
    for case in cases:
        assert set(case["sequence"]) <= set(AA)
        inp = {
            "name": case["id"],
            "dialect": "alphafold3",
            "version": 1,
            "modelSeeds": [101, 202, 303],
            "sequences": [
                {
                    "protein": {
                        "id": chain,
                        "sequence": seq,
                        "unpairedMsa": "",
                        "pairedMsa": "",
                        "templates": [],
                    }
                }
                for chain, seq in [("A", target.sequence), ("B", case["sequence"])]
            ],
        }
        path = args.output / "inputs" / f"{case['id']}.json"
        path.write_text(json.dumps(inp, indent=2) + "\n")
        case.update(
            input=str(path.relative_to(args.output)),
            input_sha256=file_hash(path),
            predictions=[],
            status="Prepared; not submitted",
            length=len(case["sequence"]),
        )
    record = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "Prepared control inputs only; no inference or submission",
        "cases": cases,
        "target": {
            "sequence": target.sequence,
            "native_range": [93, 484],
            "reference_sha256": file_hash("data/references/37NT/target.pdb"),
        },
        "protocol": {
            "engine": "AlphaFold3",
            "state": "draft",
            "model_seeds": [101, 202, 303],
            "samples_per_seed": 1,
            "num_recycles": 10,
            "templates": False,
            "msa": "single_sequence",
            "glycans_in_inference": False,
        },
        "release_requirements": [
            "At least one qualified binder nominated; evaluate it with this exact target, engine, MSA policy, seeds and sample budget.",
            "Recheck the shared GPU-hour balance; the existing AF3 reserve is not an extra control allowance.",
            "Inspect all seeds and assess native glycan context after prediction. Low confidence alone does not prove nonbinding.",
            "These JSON inputs do not use the legacy two-generator candidate-cohort launcher unchanged.",
        ],
        "sources": {
            str(p): file_hash(p)
            for p in [
                source / "O00295.json",
                source / "2H5Q_entity1.json",
                source / "CH6.json",
                Path(__file__),
            ]
        },
    }
    (args.output / "manifest.json").write_text(json.dumps(record, indent=2) + "\n")
    print("Prepared TULP2 (520 aa) and mCherry (236 aa); no predictions or jobs created")


if __name__ == "__main__":
    main()
