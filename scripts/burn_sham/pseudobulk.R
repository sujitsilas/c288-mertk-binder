# Rscript scripts/burn_sham/pseudobulk.R [input_h5ad] [output_directory]
suppressPackageStartupMessages({library(rhdf5); library(Matrix); library(limma)})
args <- commandArgs(TRUE)
default_input <- "data/burn_sham/filtered_final_burn_sham.h5ad"
if (!file.exists(default_input) && file.exists(basename(default_input))) {
  default_input <- basename(default_input)
}
f <- if(length(args)) args[1] else default_input
out <- if(length(args)>1) args[2] else "outputs/burn_sham"
if (!file.exists(f)) stop("Missing AnnData input: ", f)
dir.create(out, showWarnings=FALSE, recursive=TRUE)
annotation <- file.path(out,"mouse_peptidases_uniprot.tsv")
if (!file.exists(annotation)) {
  query <- "(organism_id:10090) AND (reviewed:true) AND (go:0008233)"
  url <- paste0("https://rest.uniprot.org/uniprotkb/stream?query=",URLencode(query,reserved=TRUE),
                "&format=tsv&fields=accession,gene_primary,protein_name,cc_subcellular_location,ec")
  temp <- tempfile(tmpdir=out)
  download.file(url,temp,mode="wb")
  stopifnot("Gene Names (primary)" %in% names(read.delim(temp,nrows=1,check.names=FALSE)))
  stopifnot(file.rename(temp,annotation))
}
save_table <- function(x,file) write.table(x,file.path(out,file),sep="\t",row.names=FALSE,quote=FALSE)
read_key <- function(key) as.vector(h5read(f,key))
category <- function(key) {
  prefix <- paste0('obs/',key,'/')
  codes <- read_key(paste0(prefix,'codes'))
  stopifnot(all(codes>=0))
  read_key(paste0(prefix,'categories'))[1L+codes]
}
meta <- data.frame(sample=category('Sample'),condition=category('Type'),day=category('Timepoint'),cell_type=category('cell_types_simple'))
genes <- read_key('var/_index')
stopifnot(!anyDuplicated(genes), nrow(unique(meta[1:3]))==length(unique(meta$sample)))
grp <- interaction(meta$sample,meta$cell_type,drop=TRUE)
units <- meta[match(levels(grp),as.character(grp)),]; units$n_cells <- tabulate(grp,nbins=nlevels(grp))
units$id <- seq_len(nrow(units))
save_table(units,'sample_cell_counts.tsv')
p <- as.numeric(read_key('X/indptr'))
n <- nrow(meta)
stopifnot(length(p)==n+1, h5readAttributes(f,'X')[['encoding-type']]=='csr_matrix')
sums <- matrix(0,length(genes),nrow(units),dimnames=list(genes,NULL))
det <- sums
lib <- numeric(n)
for (first in seq(1,n,by=768)) {
  cells <- first:min(first+767,n)
  ptr <- p[cells[1]:(tail(cells,1)+1)]
  off <- ptr[1]+1; len <- tail(ptr,1)-ptr[1]
  xx <- as.numeric(h5read(f,'X/data',start=off,count=len))
  ii <- as.integer(h5read(f,'X/indices',start=off,count=len))
  stopifnot(all(is.finite(xx)),all(xx>=0))
  m <- sparseMatrix(i=ii+1L,p=as.integer(ptr-ptr[1]),x=expm1(xx),dims=c(length(genes),length(cells)))
  lib[cells] <- colSums(m)
  g <- sparseMatrix(i=seq_along(cells),j=as.integer(grp[cells]),x=1,dims=c(length(cells),nrow(units)))
  sums <- sums+as.matrix(m%*%g)
  det <- det+as.matrix((m>0)%*%g)
}
print(summary(lib)); stopifnot(max(abs(lib-10000))<2)
ann <- read.delim(file.path(out,'mouse_peptidases_uniprot.tsv'),check.names=FALSE,quote='')
panel <- unique(na.omit(ann[['Gene Names (primary)']]))
sentinels <- c('Adam17','Mertk','Gas6','Pros1','Adam10')
pb <- function(counts,detection,d,label) {
  # Each biological sample contributes one expression profile per cell type.
  keep <- d$n_cells>=20
  d <- d[keep,,drop=FALSE]; counts <- counts[,keep,drop=FALSE]; detection <- detection[,keep,drop=FALSE]
  day_table <- table(d$day,d$condition)
  keep <- d$day %in% rownames(day_table)[apply(day_table,1,min)>=1]
  d <- d[keep,,drop=FALSE];counts<-counts[,keep,drop=FALSE];detection<-detection[,keep,drop=FALSE]
  d$condition <- factor(d$condition,levels=c('Sham','Burn'));d$day<-factor(d$day)
  stopifnot(min(table(d$condition))>=3)
  design <- model.matrix(~day+condition,d)
  expr <- sweep(counts,2,d$n_cells,'/')
  frac <- sweep(detection,2,d$n_cells,'/')
  pass <- rowSums(frac>=0.01 & detection>=3)>=3
  logexpr <- log2(expr[pass,,drop=FALSE]+0.1)
  fit <- eBayes(lmFit(logexpr,design),trend=TRUE,robust=TRUE)
  tab <- topTable(fit,coef='conditionBurn',number=Inf,sort.by='none',confint=TRUE)
  tab$gene <- rownames(tab);tab$cell_type<-label
  for (condition in levels(d$condition)){
    chosen <- d$condition==condition
    tab[[paste0('mean_',condition,'_per10k')]]<-rowMeans(expr[tab$gene,chosen,drop=FALSE])
    tab[[paste0('mean_',condition,'_detected_fraction')]]<-rowMeans(frac[tab$gene,chosen,drop=FALSE])
    tab[[paste0('n_',condition)]]<-sum(chosen)
  }
  list(result=tab, samples=d, values=expr[intersect(c(panel,sentinels),rownames(expr)),,drop=FALSE])
}
types <- sort(unique(units$cell_type))
res <- lapply(types,function(t){idx<-units$cell_type==t;pb(sums[,idx],det[,idx],units[idx,],t)})
names(res)<-types
d<-unique(meta[1:3]);d$n_cells<-as.numeric(table(factor(meta$sample,levels=d$sample)))
mix<-as.matrix(sparseMatrix(i=units$id,j=match(units$sample,d$sample),x=1))
res[['All cells']]<-pb(sums%*%mix,det%*%mix,d,'All cells')
all<-do.call(rbind,lapply(res,'[[','result'))
save_table(all,'all_genes_DE.tsv')
prot<-all[all$gene %in% panel,]
prot$panel_FDR<-NA_real_
for (pool in c(TRUE,FALSE)) {ix<-(prot$cell_type=='All cells')==pool;prot$panel_FDR[ix]<-p.adjust(prot$P.Value[ix],'BH')}
save_table(prot,'protease_DE.tsv')
save_table(all[all$gene %in% sentinels,],'receptor_ligand_sheddase.tsv')
save_table(do.call(rbind,lapply(names(res),function(t){d<-res[[t]]$samples;d$cell_type<-t;d$id<-NULL;d})),'model_samples.tsv')
per_sample<-do.call(rbind,lapply(names(res),function(t){v<-res[[t]]$values;di<-res[[t]]$samples;data.frame(gene=rep(rownames(v),ncol(v)),sample=rep(di$sample,each=nrow(v)),condition=rep(di$condition,each=nrow(v)),day=rep(di$day,each=nrow(v)),cell_type=t,mean_per10k=as.vector(v))}))
save_table(per_sample,'protease_sample_expression.tsv')
saveRDS(list(sums=sums,detected=det,groups=units,cell_totals=summary(lib)),file.path(out,'linear_pseudobulk.rds'))
writeLines(capture.output(sessionInfo()),file.path(out,'R_session.txt'))
