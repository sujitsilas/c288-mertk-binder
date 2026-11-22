"""Compare the proposed linker/slow comparator with experimental human MMP cleavage sites."""

import argparse
import csv
import itertools
import math
from collections import Counter
from pathlib import Path
from zipfile import ZipFile

AMINO = "ACDEFGHIKLMNPQRSTVWY"
CORES = ("GPQGIWGQ", "GPQGIAGQ")
PROTEASES = ("MMP1", "MMP8", "MMP9", "MMP12", "MMP14")


def profile(lines):
    sites = set()
    for line in lines.decode("ascii").splitlines():
        left, right = line.strip().split(",")
        left, right = left[-4:].rjust(4, "X"), right[:4].ljust(4, "X")
        sites.add(left + right)
    counts = [Counter(site[i] for site in sites if site[i] in AMINO) for i in range(8)]
    totals = [sum(col.values()) for col in counts]
    background = sum(counts, Counter())
    total = sum(totals)
    weights = [
        {
            aa: math.log2(((col[aa] + 0.5) / (n + 10)) / ((background[aa] + 0.5) / (total + 10)))
            for aa in AMINO
        }
        for col, n in zip(counts, totals)
    ]
    return sites, counts, weights


def analyze(archive, output):
    rows = []
    permutations = {core: sorted(set(itertools.permutations(core))) for core in CORES}
    with ZipFile(archive) as data:
        for enzyme, library in itertools.product(PROTEASES, ("G", "T")):
            key = f"WebPICS_h{enzyme}_{library}_1%/cleavage-sites.txt"
            sites, counts, weights = profile(data.read(key))

            def score(seq):
                return round(sum(w[aa] for w, aa in zip(weights, seq)), 9)

            for core in CORES:
                observed = score(core)
                null = [score(seq) for seq in permutations[core]]
                rank = (
                    100 * (sum(x < observed for x in null) + 0.5 * null.count(observed)) / len(null)
                )
                rows.append(
                    [
                        enzyme,
                        library,
                        core,
                        len(sites),
                        int(core in sites),
                        counts[5][core[5]],
                        sum(counts[5].values()),
                        f"{observed:.3f}",
                        f"{rank:.1f}",
                    ]
                )
    with output.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "human_enzyme",
                "library_G_GluC_T_trypsin",
                "P4_to_P4prime",
                "unique_windows",
                "exact_window_observed",
                "matching_P2prime_windows",
                "resolved_P2prime_windows",
                "positional_log2_enrichment",
                "within_core_permutation_percentile",
            ]
        )
        writer.writerows(rows)
    print(f"Scored {len(rows)} enzyme/library/core combinations", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    analyze(args.archive, args.output)
