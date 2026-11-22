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
vec <- function(key) as.vector(h5read(file, key))
cat_obs <- function(key) {
  prefix <- paste0("obs/", key, "/")
  vec(paste0(prefix,"categories"))[1L + vec(paste0(prefix,"codes"))]
}
save_tsv <- function(x, name) {
  write.table(x, file.path(out, name), row.names=FALSE, quote=FALSE, sep="\t")
}
genes <- vec("var/_index")
lib <- vec("obs/nCount_RNA")
intermediate_lib <- vec("obs/nCount_originalexp")
features <- vec("obs/nFeature_RNA")
sample <- factor(cat_obs("Sample"))
types <- factor(cat_obs("cell_types_simple"))
ptr <- as.double(vec("X/indptr"))
n <- length(lib)
stopifnot(!anyDuplicated(genes), length(ptr) == n+1L, all(diff(ptr) == features),
          all(lib == vec("obs/total_counts")), all(lib > 0), all(intermediate_lib > 0),
          h5readAttributes(file, "X")[["encoding-type"]] == "csr_matrix")
sample_counts <- matrix(0, length(genes), nlevels(sample), dimnames=list(genes, levels(sample)))
type_counts <- matrix(0, length(genes), nlevels(types), dimnames=list(genes, levels(types)))
detected <- numeric(length(genes))
sum_error <- integer(n)
max_round_error <- 0
max_log_error <- 0
# Compute sparse sums without allocating the full cells-by-genes matrix.
for (start in seq(1, n, by=384)) {
  cells <- start:min(start+383L, n)
  boundaries <- ptr[start:(tail(cells,1)+1L)]
  percell <- as.integer(diff(boundaries))
  offset <- boundaries[1L]
  length_chunk <- sum(percell)
  log_values <- as.double(h5read(file, "X/data", start=offset+1, count=length_chunk))
  index <- as.integer(h5read(file, "X/indices", start=offset+1, count=length_chunk))
  inner_factor <- rep(intermediate_lib[cells]/10000, times=percell)
  raw_factor <- rep(lib[cells]/10000, times=percell)
  inferred <- expm1(expm1(log_values)*inner_factor)*raw_factor
  counts <- round(inferred)
  stopifnot(all(is.finite(counts)), all(counts >= 1), all(counts <= rep(lib[cells],times=percell)))
  max_round_error <- max(max_round_error, abs(inferred-counts))
  roundtrip <- log1p(log1p(counts/raw_factor)/inner_factor)
  max_log_error <- max(max_log_error, abs(roundtrip-log_values))
  mat <- sparseMatrix(i=index+1L, p=as.integer(boundaries-offset), x=counts,
                      dims=c(length(genes),length(cells)))
  stopifnot(length(mat@x) == length_chunk)
  sum_error[cells] <- colSums(mat)-lib[cells]
  indicator <- function(group) sparseMatrix(i=seq_along(cells), j=as.integer(group[cells]),
                                              x=1, dims=c(length(cells), nlevels(group)))
  sample_counts <- sample_counts + as.matrix(mat %*% indicator(sample))
  type_counts <- type_counts + as.matrix(mat %*% indicator(types))
  detected <- detected + tabulate(index+1L, nbins=length(genes))
}
qc <- data.frame(cells=n, nonzero_entries=tail(ptr,1), total_library_counts=sum(lib),
                 cells_matching_library_totals=sum(sum_error==0),
                 cells_matching_detected_features=sum(diff(ptr)==features),
                 max_cell_total_error=max(abs(sum_error)),
                 max_pre_round_integer_distance=max_round_error,
                 max_roundtrip_logX_error=max_log_error)
print(qc)
if (any(sum_error != 0) || max_log_error > 0.0001) {
  print(which(sum_error != 0))
  stop("Reconstruction did not pass validation")
}
total <- rowSums(sample_counts)
stopifnot(sum(total)==sum(lib),all(total==rowSums(type_counts)))
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
