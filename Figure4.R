library(ggplot2)
library(cowplot)
library(gridExtra)

load("/directory/UMAP.Rdata") #computed the UMAP of each normalized matrix 
names(lUMAP)<-c("counts", "Scran", "PsiNorm", "Linnorm", "logCPM")

clust.ann<-read.csv("/directory/cluster.annotation.csv")
clust.mem<-read.csv("/directory/cluster.membership.csv")

colors<-clust.ann$cluster_color[c(8,9,59,61,63,65,66,67,70,75,28,12,15,92,29,31,93,35)]

for(i in 1:100){
  clust.mem$subclass[clust.mem$x==i]<-clust.ann$subclass_label[i]
}
clust.mem$subclass<-as.factor(clust.mem$subclass)

timesec<-c(0, 2227, 167, 258, 58)
timemin<-round(timesec/60, digits = 1)
gg2<-list()
for (i in 1:length(lUMAP)) {
  df<-data.frame(U1=lUMAP[[i]][,1],
                 U2=lUMAP[[i]][,2],
                 cell_subline=clust.mem$subclass)
  df<-df[-c(which(df$cell_subline=="doublet"),
            which(df$cell_subline=="Low Quality")), ]
  df$cell_subline<-droplevels(df$cell_subline)
  
  gg2[[i]]<-ggplot(df, aes(U1, U2, color=cell_subline))+
    geom_point(aes(color=cell_subline))+
    ggtitle(names(lUMAP)[i],
            subtitle = paste("SEC:", timesec[i], "MIN:", timemin[i]))+
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5, size = 15),
          plot.subtitle = element_text(size = 12), legend.position = "none",
          axis.title.x = element_blank(), axis.title.y = element_blank())+ 
    scale_color_manual(values = colors)
}

load("/directory/ARIPCA.Rdata") #computed the Adjusted Rand Index from PCA
rownames(ARI)[3]<-"PsiNorm"
ARI<-rbind(ARI, rep(0.527, 50))
rownames(ARI)[6]<-"ScranCLUST"
df<-data.frame(ARI=rowMeans(ARI),
               norm=rownames(ARI))
gARI<-ggplot(df, aes(x=reorder(norm, ARI), y=ARI))+
  geom_point(size=4)+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text.x = element_text(hjust = 1,size = 13, vjust = 1, angle = 30),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=15),
        axis.text.y = element_text(size = 13)
  )

plot_grid(gg2[[1]], gg2[[5]], gg2[[4]],
          gg2[[3]], gg2[[2]], gARI, nrow=2,
          labels = c("A", "B", "C", "D", "E" ,"F"),
          label_size = 10)

# scran cluster
load("/directory/UMAP2scran_clust.Rdata") #UMAP for scran with 10 core
df<-df[-c(which(df$subclass_label=="doublet"),
          which(df$subclass_label=="Low Quality")), ]
df$cell_subline<-droplevels(as.factor(df$subclass_label))

gscrac<-ggplot(df, aes(U1, U2, color=subclass_label))+
  geom_point(aes(color=subclass_label))+
  ggtitle("ScranCLUST",
          subtitle = "SEC: 2138 MIN: 35.6")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 15),
        plot.subtitle = element_text(size = 12), 
        legend.position = "none",
        axis.title.x = element_blank(), axis.title.y = element_blank())+ 
  scale_color_manual(values = colors)

legend<-get_legend(ggplot(df, aes(U1, U2, color=subclass_label))+
                     geom_point(aes(color=subclass_label))+
                     theme_classic()+
                     theme(legend.position = "bottom",
                           legend.title = element_blank())+
                     scale_color_manual(values = colors))

# HEATMAP 
sce1<-read10xCounts("/directory/umi_counts.h5")
annotation<-read.csv2(file = "/directory/cluster.annotation.csv",
                      sep = ";", header = T)
members<-read.csv2(file = "/directory/cluster.membership.csv",
                   sep = ";", header = T)
colnames(members) <- cbind("cell_id", "cluster_id")
cluster_info <- inner_join(members, annotation)# merge annotations

# split cell id in barcode + run_id
tmp <- strsplit(cluster_info$cell_id, "L")
cluster_info$barcode <- purrr::map_chr(tmp, 1)
cluster_info$run_id <- paste0("L", purrr::map_chr(tmp, 2))

# merge annotation in colData
colData(sce1) <- DataFrame(left_join(as.data.frame(colData(sce1)), 
                                     cluster_info, 
                                     by = c("Barcode" = "barcode")))
colData(sce1)
keep<-which(!is.na(sce1$cell_id))
sce1<-sce1[,keep]

scran_norm = function(sce){
  tp = system.time({
    sce = scran::computeSumFactors(sce)
    scra = scater::normalizeCounts(sce, exprs_values="counts") 
  })
  assay(sce,"Scran")<-scra
  
  method_name = "scran"
  method_type = "norm"
  if (!is.null(metadata(sce)$running_time)){
    metadata(sce)$running_time = rbind(metadata(sce)$running_time, data.frame(method=method_name, method_type=method_type, time=unname(tp)[1]))
  }else{
    metadata(sce)$running_time = data.frame(method=method_name,method_type=method_type,time=unname(tp)[1])
  }
  return(sce)
}

#PARETO
pareto.MLE2 <- function(mat) {
  n <- nrow(mat)
  m <- colMins(mat)
  a <- n/colSums(t(t(log(mat)) - log(m)))
  return(a)
}

pareto_norm2 <- function(sce) {
  c <- counts(sce)
  tp <- system.time({
    alfa <- pareto.MLE2(c+1)
    m <- log2(t(t(c)*alfa)+1)
    assay(sce, "Pareto") <- m
    sce$sizefactor<-1/alfa
  })
  method_name = "Pareto"
  method_type = "norm"
  if (!is.null(metadata(sce)$running_time)){
    metadata(sce)$running_time = rbind(metadata(sce)$running_time, data.frame(method=method_name, method_type=method_type, time=unname(tp)[1]))
  }else{
    metadata(sce)$running_time = data.frame(method=method_name,method_type=method_type,time=unname(tp)[1])
  }
  return(sce)
}

#DESeq
DESeq2_norm = function(sce){
  c<-counts(sce)
  tp = system.time({
    sizeFactors(sce) <- DESeq2::estimateSizeFactorsForMatrix(c,
                                                             type = "poscounts")
    des <-scater::normalizeCounts(sce, exprs_values="counts")
    assay(sce,"DESeq2")<-des
  })
  
  method_name = "DESeq2"
  method_type = "norm"
  if (!is.null(metadata(sce)$running_time)){
    metadata(sce)$running_time = rbind(metadata(sce)$running_time, data.frame(method=method_name, method_type=method_type, time=unname(tp)[1]))
  }else{
    metadata(sce)$running_time = data.frame(method=method_name,method_type=method_type,time=unname(tp)[1])
  }
  return(sce)
}

#TMM
TMM_norm = function(sce){
  c<-counts(sce)
  tp = system.time({
    sizeFactors(sce) <- edgeR::calcNormFactors(c, 
                                               method = "TMM") * colSums(c)
    tmmn <- scater::normalizeCounts(sce, exprs_values="counts")
    assay(sce, "TMM")<-tmmn
  })
  
  method_name = "TMM"
  method_type = "norm"
  if (!is.null(metadata(sce)$running_time)){
    metadata(sce)$running_time = rbind(metadata(sce)$running_time, data.frame(method=method_name, method_type=method_type, time=unname(tp)[1]))
  }else{
    metadata(sce)$running_time = data.frame(method=method_name,method_type=method_type,time=unname(tp)[1])
  }
  return(sce)
}

norm_method <- list(
  scran = scran_norm,
  pareto=pareto_norm2,
  DESeq2=DESeq2_norm,
  TMM=TMM_norm)

pvalb<-which(sce1$subclass_label=="Pvalb") #Oligo vs Astro
sst<-which(sce1$subclass_label=="Sst")

oligo<-which(sce1$subclass_label=="Oligo") #Oligo vs Astro
astro<-which(sce1$subclass_label=="Astro")

set.seed(1234)
pvalb<-pvalb[sample(1:length(pvalb), 400)] #subsamples
sst<-sst[sample(1:length(sst), 400)]  

counts<-cbind(counts(sce1)[,pvalb], counts(sce1)[,sst])
#counts<-cbind(counts(sce1)[,oligo], counts(sce1)[,astro])
sce<-SingleCellExperiment(assays=list(counts=as.matrix(counts)))
remove(sce1, counts)
sce$celname<-rep(c("Pvalb", "Sst"), c(length(pvalb), length(sst)))
#sce$celname<-rep(c("Oligo", "Astro"), c(length(oligo), length(astro)))
sce<-sce[rowSums(counts(sce) > 0) > 4,]

datasets<-list(dati=sce)
#normalization
res1 <- datasets %>%
  apply_methods(norm_method)

res<-res1$result
res<-SingleCellExperiment(assays=list(counts=counts(res[[1]]),
                                      Scran=assay(res[[1]], "Scran"),
                                      Pareto=assay(res[[2]], "Pareto"),
                                      DESeq2=assay(res[[3]], "DESeq2"),
                                      TMM=assay(res[[4]], "TMM")))

res$cell_type=rep(c("Pvalb", "Sst"), times=c(length(pvalb), length(sst)))
#res$cell_type=rep(c("Oligo", "Astro"), times=c(length(oligo), length(astro)))
res$cell_type=as.factor(res$cell_type)
res$SFscran<-res1$result[[1]]$sizeFactor
res$SFpareto<-res1$result[[2]]$sizefactor
res$SFdeseq<-res1$result[[3]]$sizeFactor
res$SFtmm<-res1$result[[4]]$sizeFactor
res$SFlogcpm<-1/apply(counts(res),2, sum)

library(edgeR)
names<-c(assayNames(res)[-1], "logCPM") #in 1 there is raw count matrix
sf<-colData(res)[,-1] #the size factors

qlflist<-list()
DElist<-list()
for (i in 1:length(names)) {
  y <- DGEList(counts=counts(res),
               norm.factors = sf[,i],
               group = res$cell_type)
  design <- model.matrix(~res$cell_type)
  y <- estimateDisp(y, design)
  
  fit <- glmQLFit(y, design)
  qlflist[[i]] <- glmQLFTest(fit, coef = 2)
  DElist[[i]]<-topTags(qlflist[[i]], n = nrow(qlflist[[i]]),
                       adjust.method = "BH",
                       p.value = 0.05)$table
}
names(DElist)<-names(qlflist)<-names

DE_PvalbvsSst<-DElist
for (i in 1:length(DE_PvalbvsSst)) {
  DE_PvalbvsSst[[i]]$genenames=rownames(DE_PvalbvsSst[[i]])
}

# DE_OligovsAstro<-DElist
# for (i in 1:length(DE_OligovsAstro)) {
#   DE_OligovsAstro[[i]]$genenames=rownames(DE_OligovsAstro[[i]])
# }

#once obtained the two DEG objects DE_OligovsAstro and DE_PvalbvsSst
for (i in 1:length(oligo)) {
  oli<-as_tibble(DE_OligovsAstro[[i]])
  oli<-oli %>% arrange(FDR)
  oligo[[i]]<-oli[1:nrow(oligo[[i]]),]
  
  pva<-as_tibble(DE_PvalbvsSst[[i]])
  pva<-pva %>% arrange(FDR)
  pvalb[[i]]<-pva[1:nrow(pvalb[[i]]),]
}

ensembl=useMart("ensembl")
ensembl = useDataset("mmusculus_gene_ensembl", mart=ensembl)

Glist_pvalb<-list()
Glist_oligo<-list()
for (i in 1:length(pvalb)) {
  Glist_pvalb[[i]]<-getBM(filters= "ensembl_gene_id", 
                          attributes= c("ensembl_gene_id", "mgi_symbol"),
                          values=list(pvalb[[i]]$genenames),
                          mart= ensembl)
  Glist_oligo[[i]]<-getBM(filters= "ensembl_gene_id", 
                          attributes= c("ensembl_gene_id", "mgi_symbol"),
                          values=list(oligo[[i]]$genenames),
                          mart= ensembl)
}
names(Glist_pvalb)=names(Glist_oligo)=names(pvalb)

f<-function(list1, list2){
  list3<-list1
  for (i in 1:length(list1)) {
    m<-match(list1[[i]]$genenames, list2[[i]]$ensembl_gene_id)
    list3[[i]]<-cbind(list3[[i]], list2[[i]]$mgi_symbol[m])
  }
  list3
}

DEpvalb<-f(pvalb, Glist_pvalb)
DEoligo<-f(oligo, Glist_oligo)

crearank<-function(tab,mark=markers){
  tmp<-data.frame(tab,rank=c(1:nrow(tab)))
  f<-match(mark,tab$`list2[[i]]$mgi_symbol[m]`)
  tmp2<-tmp[f,c(1,6,7)]
  tmp2
  
}

library(ComplexHeatmap)
library(cowplot)
markers<-c("Pvalb","Sst")
prova<-lapply(DEpvalb, crearank)

resPvalb<-data.frame(marker=prova[[1]]$list2..i...mgi_symbol.m.,
                     ScranFC=prova[[1]]$logFC,Scran=prova[[1]]$rank,
                     PsiNormFC=prova[[2]]$logFC,PsiNorm=prova[[2]]$rank,
                     DESEQ2FC=prova[[3]]$logFC,DESEQ2=prova[[3]]$rank,
                     TMM.FC=prova[[4]]$logFC,TMM=prova[[4]]$rank,
                     logCPM.FC=prova[[5]]$logFC,logCPM=prova[[5]]$rank,
                     ScranCLUSTFC=prova[[6]]$logFC,ScranCLUST=prova[[6]]$rank)

markers<-c("Mbp","Aqp4","Rorb")
prova<-lapply(DEoligo,crearank)
resOligo<-data.frame(marker=prova[[1]]$list2..i...mgi_symbol.m.,ScranFC=prova[[1]]$logFC,Scran=prova[[1]]$rank,
                     PsiNormFC=prova[[2]]$logFC,PsiNorm=prova[[2]]$rank,
                     DESEQ2FC=prova[[3]]$logFC,DESEQ2=prova[[3]]$rank,
                     TMM.FC=prova[[4]]$logFC,TMM=prova[[4]]$rank,
                     logCPM.FC=prova[[5]]$logFC,logCPM=prova[[5]]$rank,
                     ScranCLUSTFC=prova[[6]]$logFC,ScranCLUST=prova[[6]]$rank)

refinal<-rbind(resPvalb,resOligo)
refinal2<-refinal[,seq(3,13,by=2)]
rownames(refinal2)<-refinal[,1]
refinal2<-rbind(refinal2, "Mean"= apply(refinal2,2,mean))

library(RColorBrewer)
my_palette <- rev(colorRampPalette(brewer.pal(11,name="RdYlBu"))(50))
my_palette<-my_palette[seq(10,40,by=2)]

g1<-pheatmap(refinal2[1:2,], display_numbers = T,cluster_rows = F,cluster_cols = F,legend=F,show_colnames = F,color =my_palette, main="Pvalb vs Sst", fontsize=14)
g2<-pheatmap(refinal2[3:5,], display_numbers = T,cluster_rows = F,cluster_cols = F,legend=F,show_colnames = F,color =my_palette, main="Oligo vs Astro", fontsize=14)
g3<-pheatmap(refinal2[6,], display_numbers = T,cluster_rows = F,cluster_cols = F,color =my_palette,legend=F, main="Mean Rank", 
             fontsize = 14, angle_col = "45")

library(ggplotify)
pheat<-plot_grid(as.ggplot(g1),as.ggplot(g2),as.ggplot(g3), 
                 nrow=3, align = "v" )

## FINAL plot
gUMAP<-gridExtra::grid.arrange(gg2[[1]], gg2[[5]], gg2[[4]],
                               gg2[[3]], gg2[[2]],gscrac, ncol=3,
                               bottom="UMAP1", left="UMAP2")
gUMAP<-grid.arrange(gUMAP, bottom=legend)

pg<-plot_grid(gARI, pheat, ncol=2, 
              labels=c("B", "C"), align = "v",
              rel_widths = c(4,7))
plot_grid(gUMAP, pg, ncol=1, labels = c("A", ""), rel_heights = c(7,5))
