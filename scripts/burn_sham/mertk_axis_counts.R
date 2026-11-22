# Rscript scripts/burn_sham/mertk_axis_counts.R [input_h5ad] [output_tsv]
suppressPackageStartupMessages({library(rhdf5); library(Matrix); library(edgeR); library(limma)})

args <- commandArgs(TRUE)
file <- if (length(args)) args[1L] else "data/burn_sham/filtered_final_burn_sham.h5ad"
out <- if (length(args) > 1L) args[2L] else "outputs/burn_sham/mertk_axis_counts.tsv"
if (!file.exists(file)) stop("Missing AnnData input: ", file)
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE))
source(file.path(dirname(normalizePath(script, mustWork=TRUE)), "count_reconstruction.R"))

meta <- data.frame(sample=factor(h5_category(file, "Sample")),
                   condition=factor(h5_category(file, "Type"), levels=c("Sham", "Burn")),
                   day=factor(h5_category(file, "Timepoint")),
                   cell_type=factor(h5_category(file, "cell_types_simple")))
stopifnot(!anyNA(meta), nrow(unique(meta[1:3])) == nlevels(meta$sample))
group <- interaction(meta$sample, meta$cell_type, drop=TRUE)
profiles <- meta[match(levels(group), as.character(group)), ]
profiles$n_cells <- tabulate(group, nbins=nlevels(group))
reconstructed <- reconstruct_counts(file, list(profile=group), detect_by="profile")
stopifnot(identical(colnames(reconstructed$by_group$profile), levels(group)),
          identical(colnames(reconstructed$detected_by_group), levels(group)))

genes <- c("Mertk", "Gas6", "Pros1", "Axl", "Adam17", "Adam8", "Mmp8", "Mmp9", "Mmp12")
fit_group <- function(group_label) {
  eligible <- profiles$cell_type == group_label & profiles$n_cells >= 20
  condition_by_day <- table(profiles$day[eligible], profiles$condition[eligible])
  shared_days <- rownames(condition_by_day)[apply(condition_by_day, 1, min) >= 1]
  eligible <- eligible & profiles$day %in% shared_days
  d <- droplevels(profiles[eligible, , drop=FALSE])
  stopifnot(!anyDuplicated(d$sample), length(table(d$condition)) == 2, min(table(d$condition)) >= 3)
  design <- model.matrix(~day + condition, d)
  raw <- reconstructed$by_group$profile[, eligible, drop=FALSE]
  detected <- reconstructed$detected_by_group[, eligible, drop=FALSE]
  y <- DGEList(counts=raw)
  pass <- filterByExpr(y, design=design) &
    rowSums(detected >= 3 & sweep(detected, 2, d$n_cells, "/") >= 0.01) >= 3
  y <- calcNormFactors(y[pass, , keep.lib.sizes=FALSE], method="TMM")
  fit <- eBayes(lmFit(voom(y, design, plot=FALSE), design), robust=TRUE)
  table <- topTable(fit, coef="conditionBurn", number=Inf, sort.by="none", confint=TRUE)
  if (!all(genes %in% rownames(table))) {
    stop(group_label, " expression filter excluded: ", paste(setdiff(genes, rownames(table)), collapse=", "))
  }
  data.frame(gene=genes, cell_type=group_label, table[genes, c("logFC","CI.L","CI.R","P.Value","adj.P.Val")],
             n_Burn=sum(d$condition == "Burn"), n_Sham=sum(d$condition == "Sham"), row.names=NULL)
}

print(reconstructed$qc)
result <- do.call(rbind, lapply(c("Mφ", "Mono/MDM"), fit_group))
dir.create(dirname(out), showWarnings=FALSE, recursive=TRUE)
write.table(result, out, row.names=FALSE, quote=FALSE, sep="\t")
