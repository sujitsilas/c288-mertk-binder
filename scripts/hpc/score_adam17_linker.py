"""Apply the rounded human ADAM17 matrix in Zhou et al., Fig. 5D, to specified bonds."""

import argparse
import csv
from pathlib import Path


def run(matrix_path, receptor_path, output):
    with matrix_path.open() as handle:
        matrix = list(csv.DictReader(handle, delimiter="\t"))
    positions = ["P4", "P3", "P2", "P1", "P1prime", "P2prime", "P3prime", "P4prime"]
    if [row["position"] for row in matrix] != positions:
        raise ValueError("Matrix positions differ from the eight published subsites")
    header, *lines = receptor_path.read_text().splitlines()
    sequence = "".join(lines)
    if "|Q12866|" not in header or sequence[489:491] != "PS":
        raise ValueError("Expected the human MERTK Pro490–Ser491 bond")
    cores = {
        "TNF_reference": "LAQAVRSS",
        "APP_reference": "HHQKLVFF",
        "human_MERTK_Pro490_Ser491": sequence[486:494],
        "proposed_linker": "GPQGIWGQ",
        "alanine_comparator": "GPQGIAGQ",
    }
    scores = {
        name: round(sum(float(row[aa]) for row, aa in zip(matrix, core)), 1)
        for name, core in cores.items()
    }
    # Table 1 gives 11.56 for TNF and <0 for APP. Figure 5D prints one decimal.
    if abs(scores["TNF_reference"] - 11.56) > 0.4 or scores["APP_reference"] >= 0:
        raise ValueError("Reference scores disagree with the published control sites")
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["construct", "P4_to_P1", "P1prime_to_P4prime", "ADAM17_motif_score"])
        for name, core in cores.items():
            writer.writerow([name, core[:4], core[4:], scores[name]])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ("matrix", "receptor_fasta", "output"):
        parser.add_argument(arg, type=Path)
    args = parser.parse_args()
    run(args.matrix, args.receptor_fasta, args.output)
