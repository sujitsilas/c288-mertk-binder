# Rscript scripts/burn_sham/pooled_transcripts.R [h5ad] [output directory]
# Undo both log1p/library-normalization steps. The obs columns store their denominators.
suppressPackageStartupMessages({library(rhdf5); library(Matrix)})
args <- commandArgs(TRUE)
default_input <- "data/burn_sham/filtered_final_burn_sham.h5ad"
if (!file.exists(default_input) && file.exists(basename(default_input))) {
  default_input <- basename(default_input)
}
file <- if (length(args)) args[1] else default_input
out <- if (length(args) > 1) args[2] else "outputs/burn_sham"
if (!file.exists(file)) stop("Missing AnnData input: ", file)
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE))
source(file.path(dirname(normalizePath(script, mustWork=TRUE)), "count_reconstruction.R"))
save_tsv <- function(x, name) {
  write.table(x, file.path(out, name), row.names=FALSE, quote=FALSE, sep="\t")
}
sample <- factor(h5_category(file, "Sample"))
types <- factor(h5_category(file, "cell_types_simple"))
reconstructed <- reconstruct_counts(file, list(sample=sample, type=types))
sample_counts <- reconstructed$by_group$sample
type_counts <- reconstructed$by_group$type
detected <- reconstructed$detected
qc <- reconstructed$qc
genes <- rownames(sample_counts)
n <- qc$cells
print(qc)
total <- rowSums(sample_counts)
stopifnot(all(total == rowSums(type_counts)))
mean_sample_cpm <- rowMeans(sweep(sample_counts,2,colSums(sample_counts),"/")*1e6)
all_genes <- data.frame(gene=genes,total_transcripts=total,
                        pooled_CPM=1e6*total/sum(total),mean_sample_CPM=mean_sample_cpm,
                        cells_detected=detected,pct_cells_detected=100*detected/n,
                        samples_detected=rowSums(sample_counts>0),cell_types_detected=rowSums(type_counts>0),
                        row.names=NULL)
all_genes <- all_genes[order(-all_genes$total_transcripts, all_genes$gene),]
ann <- read.delim(file.path(out,"mouse_peptidases_uniprot.tsv"), check.names=FALSE, quote="")
proteases <- all_genes[all_genes$gene %in% ann[["Gene Names (primary)"]],]
proteases$pooled_rank <- rank(-proteases$total_transcripts, ties.method="min")
proteases$equal_sample_CPM_rank <- rank(-proteases$mean_sample_CPM, ties.method="min")
save_tsv(proteases[,c("pooled_rank",setdiff(names(proteases),"pooled_rank"))],"pooled_protease_transcripts.tsv")
save_tsv(all_genes,"pooled_gene_transcripts.tsv")
save_tsv(qc,"count_reconstruction_qc.tsv")
print(head(proteases[,c("gene","total_transcripts","pooled_rank","equal_sample_CPM_rank")],12))
print(proteases[proteases$gene %in% c("Mmp14","Mmp2","Mmp8","Mmp9","Adam17","Mmp12"),
                c("gene","total_transcripts","pooled_rank","equal_sample_CPM_rank")])
