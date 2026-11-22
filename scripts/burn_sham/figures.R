# Run after pseudobulk.R. All measured contrasts stay on a common log2 scale.
suppressPackageStartupMessages(library(ggplot2))
arg <- commandArgs(TRUE)
out <- if(length(arg)) arg[1] else "outputs/burn_sham"
tab <- read.delim(file.path(out, "protease_DE.tsv"), check.names=FALSE)
cells <- c("DC/Mono","Mono/MDM","Mφ","Neu","T","Endo.","Fibs","Krt","Seb","Smcs")
names(cells) <- c("DC / mono","Mono / MDM","Macrophage","Neutrophil","T cell","Endothelial","Fibroblast","Keratinocyte","Sebocyte","SMC")
t <- subset(tab, cell_type != "All cells")
genes <- sort(unique(t$gene))
m <- matrix(NA_real_, length(genes), length(cells), dimnames=list(genes,unname(cells)))
m[cbind(match(t$gene,genes), match(t$cell_type,cells))] <- t$logFC
# Zero is the neutral effect. Use it to assign distance only. Keep NA in every plot/table.
training <- m
training[is.na(training)] <- 0
set.seed(20261118)
k <- kmeans(training, centers=6, nstart=100, iter.max=150)
cluster_order <- order(rowMeans(k$centers))
clusters <- match(k$cluster, cluster_order)
distances <- rowSums((training - k$centers[k$cluster,,drop=FALSE])^2)
ord <- order(clusters, distances, genes)
order_genes <- genes[ord]
assignments <- data.frame(gene=order_genes, cluster=clusters[ord], tested_types=rowSums(!is.na(m))[ord], m[ord,],check.names=FALSE)
write.table(assignments,file.path(out,"heatmap_clusters.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
grid <- merge(expand.grid(gene=genes,cell_type=unname(cells),stringsAsFactors=FALSE),t,all.x=TRUE)
grid$cluster <- paste0("K", clusters[match(grid$gene,genes)])
grid$cell_type <- factor(grid$cell_type, levels=cells,labels=names(cells))
theme_set(theme_minimal(base_size=11))
heat <- function(rows, overview=FALSE, heading, subtitle) {
  d <- grid[grid$gene %in% rows,]
  d$gene <- factor(d$gene,levels=rev(rows))
  p <- ggplot(d,aes(cell_type,gene,fill=logFC)) + geom_tile() +
    geom_point(data=subset(d,!is.na(panel_FDR) & panel_FDR<.05),shape=21,
               color="#111111",fill="white",size=if(overview) 1.05 else 1.9, stroke=.45) +
    scale_fill_gradient2(low="#24588a",mid="#faf7f2",high="#b33b15",midpoint=0,
          limits=c(-2,2),oob=scales::squish,na.value="#c1c7ca",
          breaks=c(-2,-1,0,1,2),labels=c("≤−2","−1","0","+1","≥+2"),name="Burn − sham (log2)") +
    scale_y_discrete(expand=c(0,0)) + scale_x_discrete(position="top",expand=c(0,0)) +
    labs(title=heading,subtitle=subtitle,x=NULL,y=NULL,
         caption="White dot  q < 0.05 across protease × cell-type tests     |     Gray  insufficient expression for testing") +
    theme(panel.grid=element_blank(),legend.position="bottom",legend.key.width=grid::unit(1.1,"cm"),
          plot.title=element_text(face="bold",size=16),plot.subtitle=element_text(color="#394953",margin=margin(b=12)),
          axis.text.x=element_text(angle=0,hjust=.5,color="#22343c",size=8.5,margin=margin(b=5)),
          panel.spacing.y=grid::unit(2,"pt"),strip.text.y.left=element_text(angle=0,size=9),
          plot.caption=element_text(hjust=0,size=8.2),plot.margin=margin(16,18,12,16))
  if (overview) p <- p + theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
  p
}
overview <- heat(order_genes,TRUE,"Protease expression across burn and sham samples",
  sprintf("%d tested genes  •  10 cell types  •  6 k-means groups  •  day-adjusted sample-level estimates", length(genes))) +
  facet_grid(cluster~., scales="free_y",space="free_y",switch="y")
ggsave(file.path(out,"protease_heatmap.png"),overview,width=10,height=10,dpi=180,bg="white")
cairo_pdf(file.path(out,"protease_heatmap.pdf"),width=9.5,height=11,onefile=TRUE)
print(overview)
for (i in seq_along(unique(clusters))) {
  g <- order_genes[clusters[ord]==i]
  pages <- split(g,ceiling(seq_along(g)/48))
  for (j in seq_along(pages)) {
    print(heat(pages[[j]],FALSE,sprintf("Protease group K%d",i),
               sprintf("Genes %d–%d of %d  •  ordered by distance to cluster center",(j-1)*48+1,min(j*48,length(g)),length(g))))
  }
}
dev.off()
focus <- c("Mmp8","Mmp12","Adam8","Ndel1","Mmp9","Mmp13","Mmp25","Adam17","Adam10")
ggsave(file.path(out,"selected_proteases.png"),heat(focus,FALSE,"Burn-associated proteases and sheddases",
       "Selected genes  •  same color scale as the full heatmap"),width=9.5,height=4.9,dpi=180,bg="white")
samples <- read.delim(file.path(out,"protease_sample_expression.tsv"))
ad <- subset(samples,gene=="Adam17")
ordered <- c("All cells",cells)
ad$population <- factor(ad$cell_type,levels=ordered,labels=c("All cells",names(cells)))
ad$day <- factor(ad$day,levels=c("D7","D10","D14","D19"))
ad$condition <- factor(ad$condition,levels=c("Sham","Burn"))
ap <- ggplot(ad,aes(day,mean_per10k,color=condition,group=condition)) +
  geom_point(position=position_jitterdodge(seed=412,jitter.width=.18,dodge.width=.48),size=1.6,alpha=.82)+
  facet_wrap(~population,ncol=4) + scale_color_manual(values=c(Sham="#246984",Burn="#b4461f"))+
  labs(title="Adam17 expression by sample and day",subtitle="Each point represents one retained sample within a cell type. All panels share a linear scale.",
       x="Day after injury or sham procedure",y="Mean normalized expression per 10,000",color=NULL)+
  theme(legend.position="bottom",panel.grid.minor=element_blank(),strip.text=element_text(face="bold"),
        plot.title=element_text(face="bold"))
ggsave(file.path(out,"Adam17_samples.png"),ap,width=10,height=7,dpi=180,bg="white")
# Quantify abundance independently of the burn-minus-sham difference.
abundance <- t[order(t$cell_type,-t$mean_Burn_per10k),c("cell_type","gene","mean_Burn_per10k",
  "mean_Sham_per10k","mean_Burn_detected_fraction","mean_Sham_detected_fraction","logFC","panel_FDR")]
abundance$burn_expression_rank <- ave(seq_len(nrow(abundance)),abundance$cell_type,FUN=seq_along)
write.table(abundance,file.path(out,"protease_abundance_rank.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
observed <- rbind(transform(grid,condition="Burn",abundance=mean_Burn_per10k),
                  transform(grid,condition="Sham",abundance=mean_Sham_per10k))
observed$condition <- factor(observed$condition,levels=c("Sham","Burn"))
observed$gene <- factor(observed$gene,levels=rev(order_genes))
by_condition <- ggplot(observed,aes(cell_type,gene,fill=log2(1+abundance))) + geom_tile() +
  facet_grid(cluster~condition,scales="free_y",space="free_y",switch="y") +
  scale_fill_gradient(low="#fcf6e9",high="#084b56",na.value="#c1c7ca",name="log2(1 + mean per 10k)") +
  scale_y_discrete(expand=c(0,0)) + scale_x_discrete(expand=c(0,0),position="top",
     labels=c("DC/Mono","Mono/MDM","Mac","Neu","T","Endo","Fib","Krt","Seb","SMC")) +
  labs(title="Protease abundance in sham and burn",subtitle="Sample-weighted means within each cell type. Genes follow the same k-means order as the contrast heatmap.",
       x=NULL,y=NULL,caption="Gray  insufficient expression for testing   |   Cell-type and gene labels appear in the companion contrast atlas.") +
  theme(panel.grid=element_blank(),axis.text.y=element_blank(),axis.text.x=element_text(size=8),
    legend.position="bottom",plot.title=element_text(face="bold"),strip.text.y.left=element_text(angle=0))
ggsave(file.path(out,"protease_abundance.png"),by_condition,width=12,height=9,dpi=165,bg="white")
chosen <- subset(samples,cell_type %in% c("Mono/MDM","All cells") & gene %in% c("Mmp8","Mmp12","Adam17"))
leave_day <- do.call(rbind,lapply(split(chosen,list(chosen$cell_type,chosen$gene),drop=TRUE),function(x) {
  do.call(rbind,lapply(sort(unique(x$day)),function(d){
    a<-subset(x,day!=d)
    a$condition<-factor(a$condition,levels=c("Sham","Burn"))
    effect<-coef(lm(log2(mean_per10k+.1)~day+condition,a))["conditionBurn"]
    data.frame(gene=a$gene[1],cell_type=a$cell_type[1],excluded_day=d,log2_effect=unname(effect))
  }))
}))
write.table(leave_day,file.path(out,"leave_one_day_out.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
