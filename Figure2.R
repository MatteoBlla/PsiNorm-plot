library(scran)
library(scater)
library(CellBench)
library(DrImpute)
library(SAVER)
library(cluster)
library(heatmaply)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(ggrepel)
library(gridExtra)

#dataset after normalization phase
load("/directory/otherdatas_after_normalization.rdata")
norm<-c("Scran","Pareto", "DESeq2", "TMM", "logCPM", "Linnorm", "CLR", "SCT")

#transform in a list of SingleCellExperiment objects
datasets<-list()
datasets$res_cs2<-res1$result[[1]]
datasets$res_10x<-res1$result[[9]]
datasets$res_dseq<-res1$result[[17]]
datasets$res_cs51<-res1$result[[25]]
datasets$res_cs52<-res1$result[[33]]
datasets$res_cs53<-res1$result[[41]]

for (i in 2:8){
  assay(datasets$res_cs2, norm[i])<-assay(res1$result[[i]], norm[i])
  assay(datasets$res_10x, norm[i])<-assay(res1$result[[i+8]], norm[i])
  assay(datasets$res_dseq, norm[i])<-assay(res1$result[[i+16]], norm[i])
  assay(datasets$res_cs51, norm[i])<-assay(res1$result[[i+24]], norm[i])
  assay(datasets$res_cs52, norm[i])<-assay(res1$result[[i+32]], norm[i])
  assay(datasets$res_cs53, norm[i])<-assay(res1$result[[i+40]], norm[i])
}
remove(res1)
load("/directory/10x_5cl_after_normalization.rdata")
datasets$res_10x5<-res1$result[[1]]
for (i in 2:8) {
  assay(datasets$res_10x5, norm[i])<-assay(res1$result[[i]], norm [i])
}
remove(res1)

#PCA plot
pca.plot<-function(data){
  names<-assayNames(data)
  title<-c("counts","Scran","PsiNorm","DESeq2",
           "TMM","logCPM","Linnorm","CLR","sctransform")
  plotlist<-list()
  for (i in 1:9){
    data<-scater::runPCA(data, exprs_values=names[i],
                         scale=T, name=paste0("PCA.", names[i]))
    df<-data.frame(PC1=reducedDim(data, paste0("PCA.", names[i]))[,1],
                   PC2=reducedDim(data, paste0("PCA.", names[i]))[,2],
                   cell_line=data$cell_line_demuxlet
                   #cell_line=data$cluster # in base the data
    )
    plotlist[[i]]<-ggplot(df, aes(PC1, PC2, color=cell_line))+
      geom_point(aes(color=cell_line), shape=1, alpha=0.7)+ 
      stat_ellipse(type = "norm")+
      scale_color_manual(values=c("forestgreen","cyan3", "orange",
                                  "firebrick3","royalblue",
                                  "deeppink","blueviolet", "orchid",
                                  "limegreen", "cyan4","cadetblue",
                                  "brown1", "bisque4","black", 
                                  "aquamarine2", "blue4","khaki3"
      ))+
      theme_classic()+ ggtitle(title[i])+
      theme(axis.title.x = element_blank(), 
            axis.title.y = element_blank(),
            legend.position = "none",
            plot.title = element_text(size = 14))
    
  }
  gridExtra::grid.arrange(plotlist[[1]],plotlist[[5]],plotlist[[8]],
                          plotlist[[6]],plotlist[[9]],plotlist[[4]],
                          plotlist[[2]],plotlist[[3]],plotlist[[7]],
                          nrow=3, ncol=3, bottom="PC1", left="PC2"
  )
}
gpca<-pca.plot(datasets$res_cs53)

# SILHOUETTE plot
silhouette<-matrix(NA, nrow = 9, ncol = 7)
rownames(silhouette)<-c("counts", norm)
colnames(silhouette)<-names(datasets)
norm<-c("Scran","Pareto", "DESeq2", "TMM", "logCPM", "Linnorm", "CLR", "SCT")
names<-c("counts", norm)
for (i in 1:9){
  for (j in 1:7){
    data<-scater::runPCA(datasets[[j]], #compute the PCs
                         exprs_values=names[i], scale=T, #scaling each norm matrix
                         name=paste0("PCA.", names[i]))
    dist<-cluster::daisy(reducedDim(data, #compute the matrix of distances
                                    paste0("PCA.", names[i])))
    dist<-as.matrix(dist)
    silhouette[i,j]<-round(summary(
      cluster::silhouette(x=as.numeric(as.factor(data$cell_line_demuxlet)),
                          dmatrix = dist))$avg.width, digits = 3)
  }
}

load("/directory/NeMOdata.Rdata")
datasets<-Nemodata

ss<-matrix(NA, nrow = 9, ncol = 6)
colnames(ss)<-c("cell_smart", "nucleus_smart", "cell_V2", "cell_V3",
                "nucleus_V2", "nucleus_V3")
norm<-c("Scran","Pareto", "DESeq2", "TMM", "logCPM", "Linnorm", "CLR", "SCT")
names<-c("counts", norm)
for (i in 1:9){
  for (j in 1:6){
    data<-scater::runPCA(datasets[[j]], # BICCN data
                         exprs_values=names[i], scale=T,
                         name=paste0("PCA.", names[i]))
    dist<-cluster::daisy(reducedDim(data, 
                                    paste0("PCA.", names[i])))
    dist<-as.matrix(dist)
    ss[i,j]<-round(summary(
      cluster::silhouette(x=as.numeric(as.factor(data$cluster)),
                          dmatrix = dist))$avg.width, digits = 3)
  }
}

silhouette<-cbind(silhouette,ss)
remove(ss)

df<-data.frame(Silhouette=as.vector(silhouette),
               norm=rep(c("counts","Scran","PsiNorm","DESeq2",
                          "TMM","logCPM","Linnorm","CLR","sctransform"), ncol(silhouette)),
               data=c(rep("mix.CELSeq", 9),
                      rep("mix.10x", 9),
                      rep("mix.DropSeq", 9),
                      rep("mix.CELSeq51", 9),
                      rep("mix.CELSeq52", 9),
                      rep("mix.CELSeq53", 9),
                      rep("mix.10x5", 9),
                      rep("BICCN.csmart", 9),
                      rep("BICCN.nsmart", 9),
                      rep("BICCN.cV2", 9),
                      rep("BICCN.cV3", 9),
                      rep("BICCN.nV2", 9),
                      rep("BICCN.nV3", 9)),
               ncell=c(rep(" 250-300", 9),rep(" 902", 9),
                       rep(" 250-300", 9),rep(" 250-300", 9),
                       rep(" 250-300", 9),rep(" 250-300", 9),
                       rep("3918", 9),rep(" 500",54)))

df$sort<-rep(0,nrow(df))
for (i in 1:nrow(df)) {
  if(df$norm[i]=="PsiNorm"){df$sort[i]<-df$Silhouette[i]}
}

gsilh<-ggplot(df, aes(x=reorder(data, sort), y=Silhouette, group=norm))+
  geom_point(aes(color=norm,size=ncell), alpha=0.45)+  
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=15),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=15),
        legend.position = "none",
        axis.ticks.x = element_blank())+
  scale_colour_manual(values=c("#660000","#FF6600",
                               "#66FF33","#009900",
                               "#3399FF","#000099",
                               "#9900FF","#00CCFF", "#999999"))

legend<-get_legend(ggplot(df, aes(x=reorder(data, sort), y=Silhouette, group=norm))+
                     theme_classic()+
                     geom_point(aes(color=norm,size=ncell), alpha=0.45)+
                     scale_colour_manual(values=c("#660000","#FF6600",
                                                  "#66FF33","#009900",
                                                  "#3399FF","#000099",
                                                  "#9900FF","#00CCFF", "#999999"))+
                     theme(legend.text = element_text(size=12),
                           legend.title = element_blank()))

# CORRELATION plot
load("/directory/otherdatas_after_normalization.rdata")
norm<-c("Scran","Pareto", "DESeq2", "TMM", "logCPM", "Linnorm", "CLR", "SCT")

#transform in a list of SingleCellExperiment objects
datasets<-list()
datasets$res_cs2<-res1$result[[1]]
datasets$res_10x<-res1$result[[9]]
datasets$res_dseq<-res1$result[[17]]
datasets$res_cs51<-res1$result[[25]]
datasets$res_cs52<-res1$result[[33]]
datasets$res_cs53<-res1$result[[41]]

for (i in 2:8){
  assay(datasets$res_cs2, norm[i])<-assay(res1$result[[i]], norm[i])
  assay(datasets$res_10x, norm[i])<-assay(res1$result[[i+8]], norm[i])
  assay(datasets$res_dseq, norm[i])<-assay(res1$result[[i+16]], norm[i])
  assay(datasets$res_cs51, norm[i])<-assay(res1$result[[i+24]], norm[i])
  assay(datasets$res_cs52, norm[i])<-assay(res1$result[[i+32]], norm[i])
  assay(datasets$res_cs53, norm[i])<-assay(res1$result[[i+40]], norm[i])
}
remove(res1)
load("/directory/10x_5cl_after_normalization.rdata")
datasets$res_10x5<-res1$result[[1]]
for (i in 2:8) {
  assay(datasets$res_10x5, norm[i])<-assay(res1$result[[i]], norm [i])
}
remove(res1)

PC1<-list()
PC2<-list()
depth<-list()
names<-assayNames(datasets$res_cs2)
for (i in 1:length(datasets)) {
  for (j in 1:9) {
    datasets[[i]]<-scater::runPCA(datasets[[i]], # mixology data
                                  exprs_values=names[j],
                                  scale=T,
                                  name=paste0("PCA.", names[j]))
  }
}
PC1<-lapply(datasets, function(x) lapply(reducedDims(x), 
                                         function(y) y[,1]))
PC2<-lapply(datasets, function(x) lapply(reducedDims(x), 
                                         function(y) y[,2]))
depth<-lapply(datasets, function(x) apply(counts(x),2,sum))

#Correlation matrix
Corr<-matrix(NA, ncol = 9, nrow = length(datasets))
colnames(Corr)<-names
rownames(Corr)<-names(PC1)
for (i in 1:length(datasets)) {
  for (j in 1:9) {
    Corr[i,j]<-round(max(abs(cor(PC1[[i]][[j]], depth[[i]])),
                         abs(cor(PC2[[i]][[j]], depth[[i]]))), 
                     digits = 3)
  }
}

load("/directory/NeMOdata.Rdata")
datasets<-Nemodata

for (i in 1:length(datasets)) {
  for (j in 1:9) {
    datasets[[i]]<-scater::runPCA(datasets[[i]], #BICCN data
                                  exprs_values=names[j],
                                  scale=T,
                                  name=paste0("PCA.", names[j]))
  }
}
PC1<-lapply(datasets, function(x) lapply(reducedDims(x), 
                                         function(y) y[,1]))
PC2<-lapply(datasets, function(x) lapply(reducedDims(x), 
                                         function(y) y[,2]))
depth<-lapply(datasets, function(x) apply(counts(x),2,sum))

#Correlation matrix
Corr2<-matrix(NA, ncol = 9, nrow = length(datasets))
colnames(Corr2)<-names
rownames(Corr2)<-names(PC1)
for (i in 1:length(datasets)) {
  for (j in 1:9) {
    Corr2[i,j]<-round(max(abs(cor(PC1[[i]][[j]], depth[[i]])),
                          abs(cor(PC2[[i]][[j]], depth[[i]]))), 
                      digits = 3)
  }
}
Corr<-rbind(Corr,Corr2)

df<-data.frame(Correlation=as.vector(Corr),
               norm=c(rep(colnames(Corr)[1],13),
                      rep(colnames(Corr)[2],13),
                      rep("PsiNorm",13),
                      rep(colnames(Corr)[4],13),
                      rep(colnames(Corr)[5],13),
                      rep(colnames(Corr)[6],13),
                      rep(colnames(Corr)[7],13),
                      rep(colnames(Corr)[8],13),
                      rep("sctransform",13)),
               data=rep(c("mix.CELSeq", "mix.10x", "mix.DropSeq", "mix.CELSeq51",
                          "mix.CELSeq52","mix.CELSeq53","mix.10x5",
                          paste0("BICCN.",rownames(Corr)[8:13])), 9),
               ncell=rep(c(" 250-300", " 902", " 250-300"," 250-300",
                           " 250-300"," 250-300","3918", " 500",
                           " 500"," 500"," 500"," 500"," 500")))
df$sort<-rep(0,nrow(df))
for (i in 1:nrow(df)) {
  if(df$norm[i]=="PsiNorm"){df$sort[i]<-df$Correlation[i]}
}

library(ggrepel)
gcorr<-ggplot(df, aes(x=reorder(data,sort), y=Correlation, color=norm))+
  geom_point(aes(color=norm, size=ncell), alpha=0.45)+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=15),
        legend.position = "none",
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=15),
        axis.ticks.x = element_blank())+     
  scale_color_manual(values=c("#660000","#FF6600","#66FF33","#009900",
                              "#3399FF","#000099","#9900FF","#00CCFF", "#999999"))

grid.arrange(gsilh, gcorr, ncol=2)

# CONCORDANCE
load("/directory/normData.rdata") #BICCN dataset
data<-list()
data$csmart<-res[,1:500]
data$nsmart<-res[,501:1000]
data$cV2<-res[,1001:1500]
data$cV3<-res[,1501:2000]
data$nV2<-res[,2001:2500]
data$nV3<-res[,2501:3000]
remove(res, r_time)

batch=names(data)  
norm<-c("counts","Scran","Pareto", "DESeq2", "TMM", "logCPM", "Linnorm", "CLR", "SCT")
vmatrixCS=matrix(NA, 
                 nrow = nrow(data[[1]]), #matrix with nrow=numb of data, ncol=numb of norms
                 ncol = length(norm))
rownames(vmatrixCS)=rownames(data[[1]]) 
vmatrixNS=matrix(NA, 
                 nrow = nrow(data[[2]]), 
                 ncol = length(norm))
rownames(vmatrixNS)=rownames(data[[2]]) 
vmatrixc2=matrix(NA, 
                 nrow = nrow(data[[3]]), 
                 ncol = length(norm))
rownames(vmatrixc2)=rownames(data[[3]]) 
vmatrixc3=matrix(NA, 
                 nrow = nrow(data[[3]]), 
                 ncol = length(norm))
rownames(vmatrixc3)=rownames(data[[3]]) 
vmatrixn2=matrix(NA, 
                 nrow = nrow(data[[3]]), 
                 ncol = length(norm))
rownames(vmatrixn2)=rownames(data[[3]]) 
vmatrixn3=matrix(NA, 
                 nrow = nrow(data[[3]]), 
                 ncol = length(norm))
rownames(vmatrixn3)=rownames(data[[3]]) 
colnames(vmatrixCS)=colnames(vmatrixNS)=colnames(vmatrixc2)=
  colnames(vmatrixc3)=colnames(vmatrixn2)=colnames(vmatrixn3)=norm

for (j in 1:9){
  fvf1=FindVariableFeatures(assay(data[[1]], norm[j]), #compute the std gene variance
                            selection.method = "vst")
  fvf2=FindVariableFeatures(assay(data[[2]], norm[j]), 
                            selection.method = "vst")
  fvf3=FindVariableFeatures(assay(data[[3]], norm[j]), 
                            selection.method = "vst")
  fvf4=FindVariableFeatures(assay(data[[4]], norm[j]), 
                            selection.method = "vst")
  fvf5=FindVariableFeatures(assay(data[[5]], norm[j]), 
                            selection.method = "vst")
  fvf6=FindVariableFeatures(assay(data[[6]], norm[j]), 
                            selection.method = "vst")
  vmatrixCS[,j]<-fvf1[,4] #allocated in the 4th column
  vmatrixNS[,j]<-fvf2[,4]
  vmatrixc2[,j]<-fvf3[,4]
  vmatrixc3[,j]<-fvf4[,4]
  vmatrixn2[,j]<-fvf5[,4]
  vmatrixn3[,j]<-fvf6[,4]
}
m1<-matrix(NA,nrow = 500, ncol = length(norm))
m2<-matrix(NA,nrow = 500, ncol = length(norm))
m3<-matrix(NA,nrow = 500, ncol = length(norm))
m4<-matrix(NA,nrow = 500, ncol = length(norm))
m5<-matrix(NA,nrow = 500, ncol = length(norm))
m6<-matrix(NA,nrow = 500, ncol = length(norm))

for (i in 1:length(norm)) {
  m1[,i]<-names(sort(vmatrixCS[,i], decreasing = T))[1:500] #take the firsts 500 genes based on variance
  m2[,i]<-names(sort(vmatrixNS[,i], decreasing = T))[1:500]
  m3[,i]<-names(sort(vmatrixc2[,i], decreasing = T))[1:500]
  m4[,i]<-names(sort(vmatrixc3[,i], decreasing = T))[1:500]
  m5[,i]<-names(sort(vmatrixn2[,i], decreasing = T))[1:500]
  m6[,i]<-names(sort(vmatrixn3[,i], decreasing = T))[1:500]
}

conc<-matrix(NA,nrow = 15, ncol = 9)
rownames(conc)<-c("CS-NS", "CS-c2", "CS-c3", "CS-n2", "CS-n3", #between concordance 
                  "NS-c2", "NS-c3", "NS-n2", "NS-n3",
                  "c2-c3", "c2-n2", "c2-n3",
                  "c3-n2", "c3-n3",
                  "n2-n3")
colnames(conc)<-norm
for (i in 1:9) {
  conc[1,i]<-sum(m1[,i]%in%m2[,i])/length(unique(c(m1[,i], m2[,i])))# 500 # how many genes are in common between datasets
  conc[2,i]<-sum(m1[,i]%in%m3[,i])/length(unique(c(m1[,i], m3[,i])))
  conc[3,i]<-sum(m1[,i]%in%m4[,i])/length(unique(c(m1[,i], m4[,i])))
  conc[4,i]<-sum(m1[,i]%in%m5[,i])/length(unique(c(m1[,i], m5[,i])))
  conc[5,i]<-sum(m1[,i]%in%m6[,i])/length(unique(c(m1[,i], m6[,i])))
  conc[6,i]<-sum(m2[,i]%in%m3[,i])/length(unique(c(m2[,i], m3[,i])))
  conc[7,i]<-sum(m2[,i]%in%m4[,i])/length(unique(c(m2[,i], m4[,i])))
  conc[8,i]<-sum(m2[,i]%in%m5[,i])/length(unique(c(m2[,i], m5[,i])))
  conc[9,i]<-sum(m2[,i]%in%m6[,i])/length(unique(c(m2[,i], m6[,i])))
  conc[10,i]<-sum(m3[,i]%in%m4[,i])/length(unique(c(m3[,i], m4[,i])))
  conc[11,i]<-sum(m3[,i]%in%m5[,i])/length(unique(c(m3[,i], m5[,i])))
  conc[12,i]<-sum(m3[,i]%in%m6[,i])/length(unique(c(m3[,i], m6[,i])))
  conc[13,i]<-sum(m4[,i]%in%m5[,i])/length(unique(c(m4[,i], m5[,i])))
  conc[14,i]<-sum(m4[,i]%in%m6[,i])/length(unique(c(m4[,i], m6[,i])))
  conc[15,i]<-sum(m5[,i]%in%m6[,i])/length(unique(c(m1[,i], m6[,i])))
}

cmlist<-list() #within concordance
for (k in 1:10){ #split randomically ten times the data
  cm<-matrix(NA, nrow = 6, ncol=9)
  rownames(cm)<-names(data)
  colnames(cm)<-norm
  
  for (i in 1:length(data)) {
    n<-1:ncol(data[[i]])
    a<-sample(n, max(n)/2, replace = F) #split in two groups
    b<-n[-a]
    dataA<-data[[i]][,a]
    dataB<-data[[i]][,b]
    mA<-matrix(NA, nrow = 500, ncol = length(norm))
    mB<-matrix(NA, nrow = 500, ncol = length(norm))
    co<-rep(NA,9)
    names(co)<-norm
    
    for (j in 1:9){
      fvfA=FindVariableFeatures(assay(dataA, norm[j]), #std variance for each group
                                selection.method = "vst")
      fvfB=FindVariableFeatures(assay(dataB, norm[j]), 
                                selection.method = "vst")
      #matrices with variances
      vmatrixA<-fvfA[,4]
      names(vmatrixA)<-rownames(fvfA)
      vmatrixB<-fvfB[,4]
      names(vmatrixB)<-rownames(fvfB)
      
      #matrices with sorted genes by variance
      mA[,j]<-names(sort(vmatrixA, decreasing = T))[1:500] #first 500 genes with highest variance
      mB[,j]<-names(sort(vmatrixB, decreasing = T))[1:500]
      
      co[j]<-sum(mA[,j] %in% mB[,j])/length(unique(c(mA[,j], mB[,j])))#500 #within concordance
    }
    cm[i,]<-co
    
  }
  cmlist[[k]]<-cm
}

concw<-matrix(NA,nrow = 6,ncol = 9)
rownames(concw)<-rownames(cmlist[[1]])
colnames(concw)<-colnames(cmlist[[1]])
for (i in 1:6) {
  for (j in 1:9) {
    concw[i,j]<-mean(cmlist[[1]][i,j],cmlist[[2]][i,j],cmlist[[3]][i,j],
                     cmlist[[4]][i,j],cmlist[[5]][i,j],cmlist[[6]][i,j],
                     cmlist[[7]][i,j],cmlist[[8]][i,j],cmlist[[9]][i,j],
                     cmlist[[10]][i,j]) #taking the mean of the 10 replicates
  }
}

concordance<-rbind(concw, conc)

load("/directory/otherdatas_after_normalization.rdata") 
norm<-c("Scran","Pareto", "DESeq2", "TMM", "logCPM", "Linnorm", "CLR", "SCT")

data<-list() #mixology data
data$res_cs2<-res1$result[[1]]
data$res_10x<-res1$result[[9]]
data$res_dseq<-res1$result[[17]]

for (i in 2:8){
  assay(data$res_cs2, norm[i])<-assay(res1$result[[i]], norm[i])
  assay(data$res_10x, norm[i])<-assay(res1$result[[i+8]], norm [i])
  assay(data$res_dseq, norm[i])<-assay(res1$result[[i+16]], norm [i])
}
remove(res1)
batch=names(data)  
norm<-c("counts",norm)
vmatrixCS=matrix(NA, 
                 nrow = nrow(data[[1]]), 
                 ncol = length(norm))
rownames(vmatrixCS)=rownames(data[[1]]) 
vmatrix10=matrix(NA, 
                 nrow = nrow(data[[2]]), 
                 ncol = length(norm))
rownames(vmatrix10)=rownames(data[[2]]) 
vmatrixDS=matrix(NA, 
                 nrow = nrow(data[[3]]), 
                 ncol = length(norm))
rownames(vmatrixDS)=rownames(data[[3]]) 
colnames(vmatrix10)=colnames(vmatrixCS)=colnames(vmatrixDS)=norm

for (j in 1:9){
  fvfcs=FindVariableFeatures(assay(data[[1]], norm[j]), 
                             selection.method = "vst")
  fvf10=FindVariableFeatures(assay(data[[2]], norm[j]), 
                             selection.method = "vst")
  fvfds=FindVariableFeatures(assay(data[[3]], norm[j]), 
                             selection.method = "vst")
  vmatrixCS[,j]<-fvfcs[,4]
  vmatrix10[,j]<-fvf10[,4]
  vmatrixDS[,j]<-fvfds[,4]
}
cs<-matrix(NA,nrow = 500, ncol = length(norm))
sc10<-matrix(NA,nrow = 500, ncol = length(norm))
ds<-matrix(NA,nrow = 500, ncol = length(norm))

for (i in 1:length(norm)) {
  cs[,i]<-names(sort(vmatrixCS[,i], decreasing = T))[1:500]
  sc10[,i]<-names(sort(vmatrix10[,i], decreasing = T))[1:500]
  ds[,i]<-names(sort(vmatrixDS[,i], decreasing = T))[1:500]
  
}
conc<-matrix(NA,nrow = 3, ncol = 9)
rownames(conc)<-c("Cseq-10x", "Cseq-Dseq", "Dseq-10x")
colnames(conc)<-norm

for (i in 1:9) {
  conc[1,i]<-sum(cs[,i]%in%sc10[,i])/length(unique(c(cs[,i], sc10[,i])))
  conc[2,i]<-sum(cs[,i]%in%ds[,i])/length(unique(c(cs[,i], ds[,i])))
  conc[3,i]<-sum(ds[,i]%in%sc10[,i])/length(unique(c(ds[,i], sc10[,i])))
}

cmlist<-list()
for (k in 1:10) {
  cm<-matrix(NA,nrow = 3,ncol=9)
  rownames(cm)<-names(data)
  colnames(cm)<-norm
  for (i in 1:length(data)) {
    n<-1:ncol(data[[i]])
    a<-sample(n, max(n)/2, replace = F)
    b<-n[-a]
    dataA<-data[[i]][,a]
    dataB<-data[[i]][,b]
    mA<-matrix(NA, nrow = 500, ncol = length(norm))
    mB<-matrix(NA, nrow = 500, ncol = length(norm))
    co<-rep(NA,9)
    names(co)<-norm
    for (j in 1:9){
      fvfA=FindVariableFeatures(assay(dataA, norm[j]), 
                                selection.method = "vst")
      fvfB=FindVariableFeatures(assay(dataB, norm[j]), 
                                selection.method = "vst")
      
      vmatrixA<-fvfA[,4]
      names(vmatrixA)<-rownames(fvfA)
      vmatrixB<-fvfB[,4]
      names(vmatrixB)<-rownames(fvfB)
      
      mA[,j]<-names(sort(vmatrixA, decreasing = T))[1:500]
      mB[,j]<-names(sort(vmatrixB, decreasing = T))[1:500]
      
      co[j]<-sum(mA[,j]%in%mB[,j])/length(unique(c(mA[,j], mB[,j])))
    }
    cm[i,]<-co
    
  }
  cmlist[[k]]<-cm
}

concw<-matrix(NA,nrow = 3,ncol = 9)
rownames(concw)<-rownames(cmlist[[1]])
colnames(concw)<-colnames(cmlist[[1]])
for (i in 1:3) {
  for (j in 1:9) {
    concw[i,j]<-mean(cmlist[[1]][i,j],cmlist[[2]][i,j],cmlist[[3]][i,j],
                     cmlist[[4]][i,j],cmlist[[5]][i,j],cmlist[[6]][i,j],
                     cmlist[[7]][i,j],cmlist[[8]][i,j],cmlist[[9]][i,j],
                     cmlist[[10]][i,j])
  }
}

conc3cel<-rbind(concw, conc)

conc<-rbind(conc3cel[1:3,], concordance[1:6,], 
            conc3cel[4:6,], concordance[7:21,])
colnames(conc)[which(colnames(conc)=="Pareto")]="PsiNorm"
colnames(conc)[which(colnames(conc)=="SCT")]="sctransform"

library(RColorBrewer)
library(pheatmap)
rownames(conc)[1:3]<-c("CELSeq", "10x", "DropSeq")
rownames(conc)[4:9]<-c("csmart", "nsmart", "cV2", "cV3", "nV2", "nV3")
rownames(conc)[10:27]<-c("CELSeq-10x", "CELSeq-DropSeq", "DropSeq-10x",
                         "csmart-nsmart", "csmart-cV2", "csmart-cV3",
                         "csmart-nV2", "csmart-nV3", "nsmart-cV2",
                         "nsmart-cV3", "nsmart-nV2", "nsmart-nV3",
                         "cV2-cV3", "cV2-nV2", "cV2-nV3", "cV3-nV2", "cV3-nV3",
                         "nV2-nV3")
myrow<-data.frame(row.names = rownames(conc),
                  rep = c(rep("Reproducibility",9),
                          rep("Replicability", 18)),
                  data=c(rep("mixology", 3),
                         rep("BICCN", 6),
                         rep("mixology", 3),
                         rep("BICCN", 15)))

mycol=list(rep=c(Reproducibility="gold", Replicability="gold3"),
           data=c(mixology="brown1", BICCN="brown4"))

numb.color<-conc
for (i in 1:nrow(conc)) {
  for (j in 1:ncol(conc)) {
    if(conc[i,j]>0.35)  numb.color[i,j]<-"yellow2"
    else numb.color[i,j]<-"gray40"
  }}

gconc<-pheatmap(conc, cluster_rows = F, treeheight_col = 7,
                annotation_row = myrow, 
                annotation_names_row = F, fontsize_col = 12,
                annotation_colors = mycol, fontsize_row = 12,
                color=colorRampPalette(c("#F7FBFF","#FFFFD9","#C7E9B4",
                                         "#7FCDBB","#41B6C4","#1D91C0",
                                         "#225EA8","#08306B",
                                         "#081D58", "#000000"))(30),
                angle_col = 45, cutree_cols = 3,
                display_numbers = T, gaps_row = 9, legend = T,
                number_color = numb.color[,c(9,1,7,3,8,4,5,2,6)])

library(ggplotify)
pg<-plot_grid(gpca, gsilh, as.ggplot(gconc), gcorr,  ncol=2, 
              labels = c("A", "B", "", "D"))
#FIGURE 2
plot_grid(pg, legend, rel_widths = c(10,1))
