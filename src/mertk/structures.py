"""Read observed protein coordinates with explicit residue mappings."""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser
from Bio.SeqUtils import seq1

AA = "ACDEFGHIKLMNPQRSTVWY"


@dataclass
class Chain:
    sequence: str
    backbone: np.ndarray
    oxygen: np.ndarray


def read_design_chain(path, chain_id):
    """Read an explicit construct; never fill gaps or invent missing atoms."""
    parser = MMCIFParser(QUIET=True) if Path(path).suffix == ".cif" else PDBParser(QUIET=True)
    structure = parser.get_structure("design", str(path))
    if len(structure) != 1:
        raise ValueError("Choose one structural model explicitly")
    model = structure[0]
    if chain_id not in model:
        raise ValueError(f"Chain {chain_id!r} is absent from {path}")
    residues = [r for r in model[chain_id] if r.id[0] == " "]
    if not residues:
        raise ValueError("Selected chain has no protein residues")
    sequence, backbone, oxygen, mapping = [], [], [], []
    previous = None
    for residue in residues:
        aa = seq1(residue.resname, undef_code="X")
        if aa not in AA or not all(name in residue for name in ("N", "CA", "C", "O")):
            raise ValueError("Construct requires standard residues and complete N/CA/C/O backbones")
        number = residue.id[1]
        if previous is not None and number not in {previous, previous + 1}:
            raise ValueError("Numbering gaps: select a documented continuous construct first")
        previous = number
        sequence.append(aa)
        backbone.append([residue[name].coord for name in ("N", "CA", "C")])
        oxygen.append(residue["O"].coord)
        mapping.append(
            {
                "position": len(sequence),
                "author_number": number,
                "insertion_code": residue.id[2].strip(),
            }
        )
    backbone, oxygen = np.asarray(backbone), np.asarray(oxygen)
    if not np.isfinite(backbone).all() or not np.isfinite(oxygen).all():
        raise ValueError("Non-finite backbone coordinates")
    return Chain("".join(sequence), backbone, oxygen), mapping
