"""Show protease transcript totals pooled across all cells in the supplied cohort."""
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator

root = Path(__file__).resolve().parents[2]
with (root / "outputs/burn_sham/pooled_protease_transcripts.tsv").open() as source:
    measured = list(csv.DictReader(source, delimiter="\t"))
rows = [row for row in measured if int(row["total_transcripts"]) > 0]
by_gene = {row["gene"]: row for row in rows}
rows.sort(key=lambda row: int(row["pooled_rank"]))
key_genes = {"Ctsb": "#344952", "Mmp12": "#a0411f", "Adam17": "#236582", "Mmp8": "#695087"}
shown = [row["gene"] for row in rows[:3]] + ["Mmp12", "Adam17", "Mmp8"]
ink = "#20343b"

fig, (plot, table) = plt.subplots(
    1, 2, figsize=(7.25, 1.42), gridspec_kw={"width_ratios": [1, 1.6]}
)
fig.subplots_adjust(left=.095, right=.99, top=.81, bottom=.30, wspace=.29)
plot.scatter([int(row["pooled_rank"]) for row in rows],
             [int(row["total_transcripts"]) for row in rows],
             color="#a5b6bd", s=8, lw=0, rasterized=True)
for gene, color in key_genes.items():
    row = by_gene[gene]
    plot.scatter(int(row["pooled_rank"]), int(row["total_transcripts"]), s=23,
                 color=color, edgecolor="white", lw=.5, zorder=3)
plot.set(yscale="log", xlim=(-13,len(rows)+8), ylim=(20,800000))
plot.set_xticks([1,200,400,len(rows)])
plot.set_yticks([100, 10000, 500000], labels=["100", "10k", "500k"])
plot.yaxis.set_minor_locator(NullLocator())
plot.tick_params(length=2.5, colors=ink, labelsize=8, pad=2)
plot.set_xlabel("Rank by pooled count", fontsize=9, labelpad=2)
plot.set_ylabel("Transcript count", fontsize=9, labelpad=1)
plot.spines[["top","right"]].set_visible(False)
plot.spines[["bottom","left"]].set_color("#819398")
plot.set_title(f"{len(rows)} detected protease genes", loc="left", fontsize=10, weight="bold", pad=5)

table_position = table.get_position()
table.set_position([table_position.x0, .06, table_position.width, .75])
table.set(xlim=(0,1),ylim=(-.2,len(shown)+1))
table.axis("off")
table.set_title("Most abundant and design-relevant genes", loc="left", fontsize=10, weight="bold", pad=5)
positions = [.015,.36,.67,.99]
for i,(x,name) in enumerate(zip(positions,["Gene","Rank","Transcripts","Cells positive"])):
    table.text(x,len(shown)+.58,name,ha="left" if i==0 else "right",va="center",
               fontsize=8.7,weight="bold", color=ink)
table.axhline(len(shown)+.04,color="#819398",lw=.65)
for index, gene in enumerate(dict.fromkeys(shown)):
    row = by_gene[gene]
    values = [gene, row["pooled_rank"], f"{int(row['total_transcripts']):,}", f"{float(row['pct_cells_detected']):.1f}%"]
    for i,(x,value) in enumerate(zip(positions,values)):
        table.text(x,len(shown)-.51-index,value,ha="left" if i==0 else "right",va="center",
                   fontsize=8.8,color=key_genes.get(gene,ink),weight="bold" if i==0 and gene in key_genes else "normal")
fig.savefig(root / "assets/figures/burn_mask_evidence.png",dpi=320)
