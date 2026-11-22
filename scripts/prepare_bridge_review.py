#!/usr/bin/env python3
"""Prepare a sourced PS/Gla/Annexin/MERTK teaching bundle, not a docked complex."""

import argparse
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import Request, urlopen

from Bio.PDB import MMCIFParser

from mertk.provenance import file_hash

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    out = args.output.resolve()
    raw = out / "structures"
    raw.mkdir(parents=True, exist_ok=True)
    records = []
    descriptions = {
        "1NL2": "Experimental bovine prothrombin fragment + Ca + lysoPS analogue; Gla mechanism surrogate, not GAS6/PROS1",
        "1A8A": "Experimental rat Annexin V + Ca + glycerophosphoserine; headgroup analogue, no acyl tails",
        "1ANX": "Experimental human Annexin V + Ca; no PS membrane",
    }
    for pdb, description in descriptions.items():
        url = f"https://files.rcsb.org/download/{pdb}.cif"
        path = raw / f"{pdb}.cif"
        if not path.exists():
            with urlopen(
                Request(url, headers={"User-Agent": "MERTK-structure-review/1.0"}), timeout=60
            ) as r:
                path.write_bytes(r.read())
        structure = MMCIFParser(QUIET=True).get_structure(pdb, path)
        residue_names = sorted({r.resname for r in structure.get_residues()})
        records.append(
            {
                "file": path.name,
                "url": url,
                "description": description,
                "sha256": file_hash(path),
                "residue_names": residue_names,
            }
        )
    assert {"CGU", "LPS", "CA"} <= set(records[0]["residue_names"])
    assert {"GSE", "CA"} <= set(records[1]["residue_names"])

    url = "https://files.rcsb.org/ligands/download/P5S_ideal.sdf"
    lipid = raw / "P5S_ideal.sdf"
    if not lipid.exists():
        with urlopen(
            Request(url, headers={"User-Agent": "MERTK-structure-review/1.0"}), timeout=60
        ) as r:
            lipid.write_bytes(r.read())
    assert "P5S" in lipid.read_text()
    records.append(
        {
            "file": lipid.name,
            "url": url,
            "sha256": file_hash(lipid),
            "description": "CCD ideal distearoyl PS (18:0/18:0); illustrative chemical conformer, not isolated neutrophil PS or a bilayer",
        }
    )

    sources = [
        ("data/references/37NT/37NT.cif", "Experimental human MERTK, folded ectodomain only"),
        (
            "data/references/pymol_review/AF-Q14393-F1-model_v6.cif",
            "Predicted human GAS6 monomer; unmodified sequence, no native Gla chemistry or calcium",
        ),
        (
            "data/references/pymol_review/AF-P07225-F1-model_v6.cif",
            "Predicted human PROS1 monomer; unmodified sequence, no native Gla chemistry or calcium",
        ),
        (
            "data/references/pymol_review/2C5D.cif",
            "Experimental human AXL-GAS6 recognition complex; not MERTK",
        ),
    ]
    old_records = json.loads((ROOT / "data/references/pymol_review/downloads.json").read_text())
    old_urls = {r["file"]: r["url"] for r in old_records}
    for source, description in sources:
        path = ROOT / source
        shutil.copy2(path, raw / path.name)
        records.append(
            {
                "file": path.name,
                "source": source,
                "description": description,
                "url": old_urls.get(path.name, "https://files.rcsb.org/download/37NT.cif"),
                "sha256": file_hash(path),
            }
        )
    for entry in ["P07225", "Q14393"]:
        shutil.copy2(ROOT / f"data/references/pymol_review/{entry}.json", out / f"{entry}.json")
    mapping_path = ROOT / "outputs/species_review_20261111/residue_mapping.json"
    numbering = None
    if mapping_path.exists():
        row = next(r for r in json.loads(mapping_path.read_text()) if r["mouse_position"] == 485)
        assert row["human_position"] == 490 and row["human_aa"] == row["mouse_aa"] == "P"
        numbering = {
            "alignment_correspondence": row,
            "human_cleavage_bond_experimentally_established_here": False,
            "mapping_sha256": file_hash(mapping_path),
            "source": "Canonical human Q12866 / mouse Q60805 alignment prepared previously",
        }
        for entry in ["Q12866", "Q60805"]:
            source = ROOT / f"outputs/species_review_20261111/{entry}.json"
            shutil.copy2(source, out / f"species_{entry}.json")
            numbering[f"{entry}_sequence_snapshot_sha256"] = file_hash(source)
    manifest = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "purpose": "First-principles structural review; separate components and experimental partial complexes",
        "not_a_ternary_complex": True,
        "cleavage_numbering_review": numbering,
        "sources": records,
        "limitations": [
            "No neutrophil membrane was isolated, reconstructed or simulated.",
            "PS headgroup chemistry is shared across acyl variants; displayed P5S is a representative chemical model, not a neutrophil lipidomics result.",
            "1NL2 is a bovine prothrombin surrogate; its deposited CGU/Ca/lipid geometry is not a GAS6 or PROS1 pose.",
            "1A8A is rat Annexin V and a headgroup analogue, not a full lipid bilayer.",
            "Predicted GAS6/PROS1 structures do not contain post-translational gamma-carboxylation, membrane or bound Ca.",
            "Neither Annexin fusion has been structurally predicted or functionally tested here.",
        ],
    }
    (out / "sources.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Prepared {len(records)} coordinate sources in {out}")


if __name__ == "__main__":
    main()
