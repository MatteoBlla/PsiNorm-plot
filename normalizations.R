library(scater)
library(scran)
library(DESeq2)
library(scran)
library(edgeR)
library(Linnorm)
library(Seurat)
library(magrittr)
library(CellBench)
NUM_OF_THREAD=8

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
    sizeFactors(sce) <- edgeR::calcNormFactors(c, method = "TMM")
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
