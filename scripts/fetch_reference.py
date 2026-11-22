#!/usr/bin/env python3
"""Fetch experimental human MERTK coordinates and export the observed construct."""

import argparse
import json
from pathlib import Path

import requests
from Bio.PDB import PDBIO, MMCIFParser, Select

from mertk.provenance import file_hash
from mertk.structures import read_design_chain


class ProteinChain(Select):
    def accept_chain(self, chain):
        return chain.id == "A"

    def accept_residue(self, residue):
        return residue.id[0] == " "


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("data/references/37NT"))
    args = parser.parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        raise FileExistsError("Choose an empty reference directory")
    url = "https://files.rcsb.org/download/37NT.cif"
    response = requests.get(url, timeout=60)
    response.raise_for_status()
    args.output.mkdir(parents=True, exist_ok=True)
    source = args.output / "37NT.cif"
    source.write_bytes(response.content)
    chain, mapping = read_design_chain(source, "A")
    structure = MMCIFParser(QUIET=True).get_structure("37NT", str(source))
    writer = PDBIO()
    writer.set_structure(structure)
    target = args.output / "target.pdb"
    writer.save(str(target), select=ProteinChain())
    manifest = {
        "pdb_id": "37NT",
        "uniprot": "Q12866",
        "source_url": url,
        "source_sha256": file_hash(source),
        "target_sha256": file_hash(target),
        "author_chain": "A",
        "sequence": chain.sequence,
        "residue_mapping": mapping,
        "construct": "Observed residues only; no missing residues reconstructed",
        "glycans": "Retained in source mmCIF; omitted from protein-only target PDB",
        "required_review": [
            "glycan accessibility",
            "ligand overlap",
            "membrane placement",
            "species-specific cleavage mapping",
            "epitope choice",
        ],
        "epitope_selected": False,
    }
    (args.output / "provenance.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Saved {len(chain.sequence)} observed residues and mapping to {args.output}")
    print("Review glycans and epitope before design; target.pdb contains protein only.")


if __name__ == "__main__":
    main()
