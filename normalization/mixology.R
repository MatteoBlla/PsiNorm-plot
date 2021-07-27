# Normalization Tian datasets

library(DESeq2)
library(scran)
library(edgeR)
library(Linnorm)
library(Seurat)
library(CellBench)
library(DrImpute)
library(SAVER)
NUM_OF_THREAD=8

load("/directory/sincell_with_class.rdata")
load("/directory/sincell_with_class_5cl.rdata")

datasets <- list(
  sc_CELseq2=sce_sc_CELseq2_qc,
  sc_10x=sce_sc_10x_qc,
  sc_Droseq=sce_sc_Dropseq_qc,
  sc_Celseq2_5cl_p1=sc_Celseq2_5cl_p1,
  sc_Celseq2_5cl_p2=sc_Celseq2_5cl_p2,
  sc_Celseq2_5cl_p3=sc_Celseq2_5cl_p3
  #sc_10x_5cl=sce_sc_10x_5cl_qc #too large I work on it separately
)
remove(sce_sc_CELseq2_qc,
       sce_sc_10x_qc,
       sce_sc_Dropseq_qc,
       sce_sc_10x_5cl_qc,
       sc_Celseq2_5cl_p1,
       sc_Celseq2_5cl_p2,
       sc_Celseq2_5cl_p3
       #sce_sc_10x_5cl_qc
       )

#filtering
geni_10x<-rowSums(counts(datasets$sc_10x)>0)>4
geni_cseq2<-rowSums(counts(datasets$sc_CELseq2)>0)>4
geni_dseq<-rowSums(counts(datasets$sc_Droseq)>0)>4 #unnecessary in 5cl datasets

datasets$sc_10x<-datasets$sc_10x[geni_10x,]
datasets$sc_CELseq2<-datasets$sc_CELseq2[geni_cseq2,]
datasets$sc_Droseq<-datasets$sc_Droseq[geni_dseq,]
message(paste0("geni 10X filtrati: ", sum(geni_10x==0),
               "// geni CELSeq2 filtrati: ", sum(geni_cseq2==0),
               "// geni DropSeq filtrati: ", sum(geni_dseq==0)))

remove(geni_10x, geni_cseq2,geni_dseq)

#normalization methods
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
library(DelayedArray)
library(matrixStats)

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

#CLR
clr_norm<-function(sce){
  c<-counts(sce)
  clr<-c*0
  tp=system.time({
    clr<-Seurat::NormalizeData(c, 
                               margin=2, normalization.method="CLR")
  })
  assay(sce, "CLR")<-clr  
  method_name = "CLR"
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
    sizeFactors(sce) <- DESeq2::estimateSizeFactorsForMatrix(c)
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
    sizeFactors(sce) <- edgeR::calcNormFactors(c, method = "TMM")* colSums(counts(sce))
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

#logCPM
logCPM_norm = function(sce){
  c<-counts(sce)
  tp = system.time({
    assay(sce, "logCPM") = log2(edgeR::cpm(c) + 1)
  })
  
  method_name = "logCPM"
  method_type = "norm"
  if (!is.null(metadata(sce)$running_time)){
    metadata(sce)$running_time = rbind(metadata(sce)$running_time, data.frame(method=method_name, method_type=method_type, time=unname(tp)[1]))
  }else{
    metadata(sce)$running_time = data.frame(method=method_name,method_type=method_type,time=unname(tp)[1])
  }
  return(sce)
}

#LINNORM
linnorm_norm = function(sce){
  c<-counts(sce)
  tp = system.time({
    assay(sce, "Linnorm") = Linnorm(c)
  })
  
  method_name = "Linnorm"
  method_type = "norm"
  if (!is.null(metadata(sce)$running_time)){
    metadata(sce)$running_time = rbind(metadata(sce)$running_time, data.frame(method=method_name, method_type=method_type, time=unname(tp)[1]))
  }else{
    metadata(sce)$running_time = data.frame(method=method_name,method_type=method_type,time=unname(tp)[1])
  }
  return(sce)
}

#SCT
SCT_norm<-function(sce){
  c<-counts(sce)
  tp = system.time({
    pbmc <- CreateSeuratObject(counts = c)
    pbmc <- SCTransform(object = pbmc, verbose = FALSE)
    assay(sce,"SCT")<-as.matrix(pbmc@assays$SCT@data)
  })
  method_name = "SCT"
  method_type = "norm"
  if (!is.null(metadata(sce)$running_time)){
    metadata(sce)$running_time = rbind(metadata(sce)$running_time,
                                       data.frame(method=method_name, method_type=method_type, time=unname(tp)[1]))
  }else{
    metadata(sce)$running_time = data.frame(method=method_name,
                                            method_type=method_type, time=unname(tp)[1])
  }
  return(sce)
}

norm_method <- list(
  scran = scran_norm,
  pareto=pareto_norm2,
  DESeq2=DESeq2_norm,
  TMM=TMM_norm,
  logCPM=logCPM_norm,
  Linnorm=linnorm_norm,
  CLR=clr_norm,
  SCT=SCT_norm
)

#normalization
res1 <- datasets %>%
  apply_methods(norm_method)
#check
a<-NULL
for (i in 1:nrow(res1)) {
  a[i]<-class(res1$result[[i]])=="SingleCellExperiment"
}
if(sum(a)==nrow(res1)){print("Normalizations have been successful for each dataset")} #if is not equal to nrow(res1) there is a problem
remove(a)

save(res1, file="/directory/otherdatas_after_normalization.Rdata")
#save(res1, file="/directory/10x_5cl_after_normalization.Rdata")


