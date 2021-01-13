library(scran)
library(scater)
library(CellBench)
library(DrImpute)
library(SAVER)
library(cluster)
library(heatmaply)
library(ggplot2)
library(RColorBrewer)

#dataset after normalization phase
load("C:/Users/Matteo/Desktop/Borsa studio/data_norms/otherdatas_after_normalization.rdata")
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
  assay(datasets$res_10x, norm[i])<-assay(res1$result[[i+8]], norm [i])
  assay(datasets$res_dseq, norm[i])<-assay(res1$result[[i+16]], norm [i])
  assay(datasets$res_cs51, norm[i])<-assay(res1$result[[i+24]], norm [i])
  assay(datasets$res_cs52, norm[i])<-assay(res1$result[[i+32]], norm [i])
  assay(datasets$res_cs53, norm[i])<-assay(res1$result[[i+40]], norm [i])
}
remove(res1)
load("C:/Users/Matteo/Desktop/Borsa studio/data_norms/10x_5cl_after_normalization.rdata")
datasets$res_10x5<-res1$result[[1]]
for (i in 2:8) {
  assay(datasets$res_10x5, norm[i])<-assay(res1$result[[i]], norm [i])
}
remove(res1)

#PCA plot
pca.plot<-function(data){
  names<-assayNames(data)
  names<-names[-2]
  title<-c("counts","Scran","PsiNorm","DESeq2","TMM","logCPM","Linnorm","CLR","sctransform")
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
            plot.title = element_text(size = 10))
    
  }
  gridExtra::grid.arrange(plotlist[[1]],plotlist[[5]],plotlist[[8]],
                          plotlist[[6]],plotlist[[9]],plotlist[[4]],
                          plotlist[[2]],plotlist[[3]],plotlist[[7]],
                          nrow=3, ncol=3,
                          bottom="PC1", left="PC2")
}
gpca<-pca.plot(datasets$res_cs53)

#create the silhouette matrix
silhouette<-matrix(NA, nrow = 9, ncol = 7)
rownames(silhouette)<-c("counts", norm)
colnames(silhouette)<-names(datasets)
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
      silhouette(x=as.numeric(as.factor(data$cell_line_demuxlet)),
                 dmatrix = dist))$avg.width, digits = 3)
  }
}

#NeMo dataset
load("C:/Users/Matteo/Desktop/Fedez data/RData/normData.rdata")
datasets<-list()
datasets$sce1<-res[,1:500]
datasets$sce2<-res[,501:1000]
datasets$sce3<-res[,1001:1500]
datasets$sce4<-res[,1501:2000]
datasets$sce5<-res[,2001:2500]
datasets$sce6<-res[,2501:3000]
remove(res, r_time)

ss<-matrix(NA,nrow = 9, ncol = 6)
colnames(ss)<-c("cell_smart", "nucleus_smart", "cell_V2", "cell_V3",
                "nucleus_V2", "nucleus_V3")
for (i in 1:9){
  for (j in 1:6){
    data<-scater::runPCA(datasets[[j]],
                         exprs_values=names[i], scale=T,
                         name=paste0("PCA.", names[i]))
    dist<-cluster::daisy(reducedDim(data, 
                                    paste0("PCA.", names[i])))
    dist<-as.matrix(dist)
    ss[i,j]<-round(summary(
      silhouette(x=as.numeric(as.factor(data$cluster)),
                 dmatrix = dist))$avg.width, digits = 3)
  }
}

silhouette<-cbind(silhouette,ss)
remove(ss)

df<-data.frame(Silhouette=as.vector(silhouette),
               norm=rep(c("counts","Scran","PsiNorm","DESeq2",
                          "TMM","logCPM","Linnorm","CLR","sctransform"), ncol(silhouette)),
               data=c(rep("CELSeq", 9),
                      rep("10x", 9),
                      rep("DropSeq", 9),
                      rep("CELSeq51", 9),
                      rep("CELSeq52", 9),
                      rep("CELSeq53", 9),
                      rep("10x5", 9),
                      rep("csmart", 9),
                      rep("nsmart", 9),
                      rep("cV2", 9),
                      rep("cV3", 9),
                      rep("nV2", 9),
                      rep("nV3", 9)),
               ncell=c(rep(" 250-300", 9),rep(" 902", 9),
                       rep(" 250-300", 9),rep(" 250-300", 9),
                       rep(" 250-300", 9),rep(" 250-300", 9),
                       rep("3918", 9),rep(" 500",54)))
df$nu[df$norm=="PsiNorm"]<-"PsiNorm"
df$nu[df$norm=="Scran"]<-"Scran"
df$nu[df$norm=="logCPM"]<-"logCPM"
df$sort<-rep(0,nrow(df))
for (i in 1:nrow(df)) {
  if(df$norm[i]=="PsiNorm"){df$sort[i]<-df$Silhouette[i]}
}

library(ggplot2)
library(ggrepel)
g1<-ggplot(df, aes(x=reorder(data, sort), y=Silhouette, group=norm))+
  #geom_text(aes(label = nu, color=norm), size=3)+
  geom_text_repel(aes(label=nu, color=norm), size=2)+
  geom_point(aes(color=norm,size=ncell), alpha=0.45)+  
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=13),
        legend.position = "none",
        axis.ticks.x = element_blank())+
  scale_colour_manual(values=c("#660000","#FF6600",
                               "#66FF33","#009900",
                               "#3399FF","#000099",
                               "#9900FF","#00CCFF", "#999999"))
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend<-get_legend(ggplot(df, aes(x=reorder(data, sort), y=Silhouette, group=norm))+
                     theme_classic()+
                     geom_point(aes(color=norm,size=ncell), alpha=0.45)+
                     scale_colour_manual(values=c("#660000","#FF6600",
                                                  "#66FF33","#009900",
                                                  "#3399FF","#000099",
                                                  "#9900FF","#00CCFF", "#999999"))+
                     theme(legend.text = element_text(size=8))
)
g.silh<-gridExtra::grid.arrange(g1, nrow=1, ncol=1, right=legend)

gridExtra::grid.arrange(gpca, g.silh, ncol=2, widths=c(0.45,0.55)) #Figure 3

#for supplementary
pca.plot(datasets$res_cs2)
pca.plot(datasets$res_10x)
pca.plot(datasets$res_dseq)
pca.plot(datasets$res_cs51)
pca.plot(datasets$res_cs52)
pca.plot(datasets$res_10x5)

