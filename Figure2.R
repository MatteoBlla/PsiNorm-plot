library(scran)
library(edgeR)
library(ggplot2)
library(magrittr)
library(CellBench)
library(scater)
library(gridExtra)

#loading the datasets after the normalizations
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

# counts and pareto normalized matrices
datas<-lapply(datasets, function(x) counts(x)) 
datas_pareto<-lapply(datasets, function(x) assay(x, "Pareto")) 

#computing slope and intercepts of the linear model between ranks and expressions
lm.reads.ranghi<-function(data){
  param<-list()
  for (j in 1:length(data)){
    slop<-c()
    inter<-c()
    for (i in 1:ncol(data[[j]])){
      rank1<-rank(unique((sort(data[[j]][,i], decreasing = T))))
      u.values1<-unique(sort(data[[j]][,i]))
      sl<-summary(lm(log1p(u.values1)~log1p(rank1)))$coefficients[2] #slopes
      int<-summary(lm(log1p(u.values1)~log1p(rank1)))$coefficients[1] #intercepts
      
      slop<-c(slop,sl)
      inter<-c(inter,int)
    }
    param$slope<-c(param$slope, slop)
    param$intercept<-c(param$intercept, inter)
  }
  param
}

param<-lm.reads.ranghi(datas)
param.pareto<-lm.reads.ranghi(datas_pareto)

#preparing the plot
# df<-data.frame(slope=c(param$slope, param.pareto$slope),
#                intercept=c(param$intercept, param.pareto$intercept),
#                data=rep(c(rep(names(datas)[1], ncol(datas[[1]])),
#                           rep(names(datas)[2], ncol(datas[[2]])),
#                           rep(names(datas)[3], ncol(datas[[3]])),
#                           rep(names(datas)[4], ncol(datas[[4]])),
#                           rep(names(datas)[5], ncol(datas[[5]])),
#                           rep(names(datas)[6], ncol(datas[[6]])),
#                           rep(names(datas)[7], ncol(datas[[7]]))),2),
#                raw=c(rep("raw", length(param$slope)),
#                      rep("normPareto", length(param$slope))))

param$data<-c(rep("cseq2", ncol(datas$res_cs2)),
              rep("10x", ncol(datas$res_10x)),
              rep("dseq", ncol(datas$res_dseq)),
              rep("cseq51", ncol(datas$res_cs51)),
              rep("cseq52", ncol(datas$res_cs52)),
              rep("cseq53", ncol(datas$res_cs53)),
              rep("10x5", ncol(datas$res_10x5))
)
slope.list<-list()
int.list<-list()

for (i in 1:7) {
  df<-data.frame(slope=c(param$slope[param$data==unique(param$data)[i]],
                         param.pareto$slope[param$data==unique(param$data)[i]]),               intercept=c(param$intercept[param$data==unique(param$data)[i]],                           param.pareto$intercept[param$data==unique(param$data)[i]]),
                 data=c(rep("Raw", ncol(datas[[i]])),
                        rep("PsiNorm", ncol(datas[[i]]))))
  p<-ggplot(df)+
    geom_histogram(aes(x=slope, fill=data),
                   position = "identity", alpha=0.6)+
    theme_classic()+
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank())+
    scale_x_continuous(breaks = c(-1.5, -1, -0.5))
  slope.list[[i]]<-p
  q<-ggplot(df)+
    geom_histogram(aes(x=intercept, fill=data),
                   position = "identity", alpha=0.6)+
    theme_classic()+
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  int.list[[i]]<-q
}

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend<-get_legend(ggplot(df)+
                     geom_histogram(aes(x=intercept, fill=data),
                                    position = "identity", alpha=0.6)+
                     theme_minimal(base_line_size = .75)+ 
                     theme(legend.position = "right",
                           legend.title = element_blank()))
library(gridExtra)
t1<-textGrob("CelSeq")
t2<-textGrob("10x")
t3<-textGrob("DropSeq")
t4<-textGrob("CelSeq51")
t5<-textGrob("CelSeq52")
t6<-textGrob("CelSeq53")
t7<-textGrob("10x5")
r1<-textGrob("Slope")
r2<-textGrob("Intercept")

graf1<-grid.arrange(right=legend, 
                    r1, slope.list[[1]],slope.list[[2]],slope.list[[3]],
                    r2, int.list[[1]],int.list[[2]],int.list[[3]],
                    ncol=4, nrow=2, 
                    bottom="Estimation", left="Absolute Frequencies")
graf2<-grid.arrange(right=legend, 
                    r1, slope.list[[4]],slope.list[[5]],slope.list[[6]],
                    slope.list[[7]],
                    r2, int.list[[4]],
                    int.list[[5]],int.list[[6]],int.list[[7]],
                    ncol=5, nrow=2, 
                    bottom="Estimation", left="Absolute Frequencies")

# taking the cells with min, median and max sequencing depth
depth<-lapply(datas, function(x) apply(x,2,sum))
depth<-lapply(depth, sort)
cell.depth<-lapply(depth, function(x) 
  c(names(x[1]), names(x[round(length(x)/2)]), names(x[length(x)])))

ff<-function(list1, list2){
  list3<-list()
  for (i in 1:length(list1)) {
    list3[[i]]<-list1[[i]][,list2[[i]]]
  }
  names(list3)<-names(list1)
  list3
}

counts_3cel<-ff(datas, cell.depth) 
par_3cel<-ff(datas_pareto, cell.depth)

#------------ would need a refactoring ------------
#CSEQ
counts<-counts_3cel$res_cs2
rank<-c()
u.val<-c()
for (i in 1:ncol(counts)) {
  rank1<-rank(unique((sort(counts[,i], decreasing = T))))
  u.values1<-unique(sort(counts[,i]))
  rank<-c(rank,rank1)
  u.val<-c(u.val,u.values1)
}

df<-data.frame(rank=rank, u.values=u.val, 
               cell_expression=c(rep("low",
                                     length((unique((sort(counts[,1])))))),
                                 rep("med", length((unique((sort(counts[,2])))))),
                                 rep("high", length((unique((sort(counts[,3]))))))
               ))
g1<-ggplot(df, aes(x=log1p(rank), y=log1p(u.values),
                   color=cell_expression))+
  geom_point(shape=1)+
  geom_smooth(method="lm", fullrange=T, se=FALSE)+ 
  scale_color_manual(values=c('#56B4E9','#E69F00', '#999999'))+
  theme_minimal(base_line_size = .75)+ 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold",hjust = 0.5))+
  ggtitle("CELSeq")

pareto<-par_3cel$res_cs2
rank<-c()
u.val<-c()
for (i in 1:ncol(pareto)) {
  rank1<-rank(unique((sort(pareto[,i], decreasing = T))))
  u.values1<-unique(sort(pareto[,i]))
  rank<-c(rank,rank1)
  u.val<-c(u.val,u.values1)
}
df<-data.frame(rank=rank, u.values=u.val, 
               cell_expression=c(rep("low",
                                     length((unique((sort(pareto[,1])))))),
                                 rep("med", length((unique((sort(pareto[,2])))))),
                                 rep("high", length((unique((sort(pareto[,3]))))))
               ))
g2<-ggplot(df, aes(x=log1p(rank), y=log1p(u.values),
                   color=cell_expression))+
  geom_point(shape=1)+
  geom_smooth(method="lm", fullrange=T, se=FALSE)+ 
  scale_color_manual(values=c('#56B4E9','#E69F00', '#999999'))+
  theme_minimal(base_line_size = .75)+ 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())


#10x
counts<-counts_3cel$res_10x
rank<-c()
u.val<-c()
for (i in 1:ncol(counts)) {
  rank1<-rank(unique((sort(counts[,i], decreasing = T))))
  u.values1<-unique(sort(counts[,i]))
  rank<-c(rank,rank1)
  u.val<-c(u.val,u.values1)
}
df<-data.frame(rank=rank, u.values=u.val, 
               cell_expression=c(rep("low",
                                     length((unique((sort(counts[,1])))))),
                                 rep("med", length((unique((sort(counts[,2])))))),
                                 rep("high", length((unique((sort(counts[,3]))))))
               ))
g5<-ggplot(df, aes(x=log1p(rank), y=log1p(u.values),
                   color=cell_expression))+
  geom_point(shape=1)+
  geom_smooth(method="lm", fullrange=T, se=FALSE)+ 
  scale_color_manual(values=c('#56B4E9','#E69F00', '#999999'))+
  theme_minimal(base_line_size = .75)+ 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold",hjust = 0.5))+
  ggtitle("10x")

pareto<-par_3cel$res_10x
rank<-c()
u.val<-c()
for (i in 1:ncol(pareto)) {
  rank1<-rank(unique((sort(pareto[,i], decreasing = T))))
  u.values1<-unique(sort(pareto[,i]))
  rank<-c(rank,rank1)
  u.val<-c(u.val,u.values1)
}
df<-data.frame(rank=rank, u.values=u.val, 
               cell_expression=c(rep("low",
                                     length((unique((sort(pareto[,1])))))),
                                 rep("med", length((unique((sort(pareto[,2])))))),
                                 rep("high", length((unique((sort(pareto[,3]))))))
               ))
g6<-ggplot(df, aes(x=log1p(rank), y=log1p(u.values),
                   color=cell_expression))+
  geom_point(shape=1)+
  geom_smooth(method="lm", fullrange=T, se=FALSE)+ 
  scale_color_manual(values=c('#56B4E9','#E69F00', '#999999'))+
  theme_minimal(base_line_size = .75)+ 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

#dseq
counts<-counts_3cel$res_dseq
rank<-c()
u.val<-c()
for (i in 1:ncol(counts)) {
  rank1<-rank(unique((sort(counts[,i], decreasing = T))))
  u.values1<-unique(sort(counts[,i]))
  rank<-c(rank,rank1)
  u.val<-c(u.val,u.values1)
}
df<-data.frame(rank=rank, u.values=u.val, 
               cell_expression=c(rep("low",
                                     length((unique((sort(counts[,1])))))),
                                 rep("med", length((unique((sort(counts[,2])))))),
                                 rep("high", length((unique((sort(counts[,3]))))))
               ))
g7<-ggplot(df, aes(x=log1p(rank), y=log1p(u.values),
                   color=cell_expression))+
  geom_point(shape=1)+
  geom_smooth(method="lm", fullrange=T, se=FALSE)+ 
  scale_color_manual(values=c('#56B4E9','#E69F00', '#999999'))+
  theme_minimal(base_line_size = .75)+ 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold",hjust = 0.5))+
  ggtitle("DropSeq")

pareto<-par_3cel$res_dseq
rank<-c()
u.val<-c()
for (i in 1:ncol(pareto)) {
  rank1<-rank(unique((sort(pareto[,i], decreasing = T))))
  u.values1<-unique(sort(pareto[,i]))
  rank<-c(rank,rank1)
  u.val<-c(u.val,u.values1)
}
df<-data.frame(rank=rank, u.values=u.val, 
               cell_expression=c(rep("low",
                                     length((unique((sort(pareto[,1])))))),
                                 rep("med", length((unique((sort(pareto[,2])))))),
                                 rep("high", length((unique((sort(pareto[,3]))))))
               ))
g8<-ggplot(df, aes(x=log1p(rank), y=log1p(u.values),
                   color=cell_expression))+
  geom_point(shape=1)+
  geom_smooth(method="lm", fullrange=T, se=FALSE)+ 
  scale_color_manual(values=c('#56B4E9','#E69F00', '#999999'))+
  theme_minimal(base_line_size = .75)+ 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

#CSEQ1
counts<-counts_3cel$res_cs51
rank<-c()
u.val<-c()
for (i in 1:ncol(counts)) {
  rank1<-rank(unique((sort(counts[,i], decreasing = T))))
  u.values1<-unique(sort(counts[,i]))
  rank<-c(rank,rank1)
  u.val<-c(u.val,u.values1)
}
df<-data.frame(rank=rank, u.values=u.val, 
               cell_expression=c(rep("low",
                                     length((unique((sort(counts[,1])))))),
                                 rep("med", length((unique((sort(counts[,2])))))),
                                 rep("high", length((unique((sort(counts[,3]))))))
               ))
g9<-ggplot(df, aes(x=log1p(rank), y=log1p(u.values),
                   color=cell_expression))+
  geom_point(shape=1)+
  geom_smooth(method="lm", fullrange=T, se=FALSE)+ 
  scale_color_manual(values=c('#56B4E9','#E69F00', '#999999'))+
  theme_minimal(base_line_size = .75)+ 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold",hjust = 0.5))+
  ggtitle("CELSeq51")

pareto<-par_3cel$res_cs51
rank<-c()
u.val<-c()
for (i in 1:ncol(pareto)) {
  rank1<-rank(unique((sort(pareto[,i], decreasing = T))))
  u.values1<-unique(sort(pareto[,i]))
  rank<-c(rank,rank1)
  u.val<-c(u.val,u.values1)
}
df<-data.frame(rank=rank, u.values=u.val, 
               cell_expression=c(rep("low",
                                     length((unique((sort(pareto[,1])))))),
                                 rep("med", length((unique((sort(pareto[,2])))))),
                                 rep("high", length((unique((sort(pareto[,3]))))))
               ))
g10<-ggplot(df, aes(x=log1p(rank), y=log1p(u.values),
                    color=cell_expression))+
  geom_point(shape=1)+
  geom_smooth(method="lm", fullrange=T, se=FALSE)+ 
  scale_color_manual(values=c('#56B4E9','#E69F00', '#999999'))+
  theme_minimal(base_line_size = .75)+ 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

#CSEQ52
counts<-counts_3cel$res_cs52
rank<-c()
u.val<-c()
for (i in 1:ncol(counts)) {
  rank1<-rank(unique((sort(counts[,i], decreasing = T))))
  u.values1<-unique(sort(counts[,i]))
  rank<-c(rank,rank1)
  u.val<-c(u.val,u.values1)
}
df<-data.frame(rank=rank, u.values=u.val, 
               cell_expression=c(rep("low",
                                     length((unique((sort(counts[,1])))))),
                                 rep("med", length((unique((sort(counts[,2])))))),
                                 rep("high", length((unique((sort(counts[,3]))))))
               ))
g11<-ggplot(df, aes(x=log1p(rank), y=log1p(u.values),
                    color=cell_expression))+
  geom_point(shape=1)+
  geom_smooth(method="lm", fullrange=T, se=FALSE)+ 
  scale_color_manual(values=c('#56B4E9','#E69F00', '#999999'))+
  theme_minimal(base_line_size = .75)+ 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold",hjust = 0.5))+
  ggtitle("CELSeq52")

pareto<-par_3cel$res_cs52
rank<-c()
u.val<-c()
for (i in 1:ncol(pareto)) {
  rank1<-rank(unique((sort(pareto[,i], decreasing = T))))
  u.values1<-unique(sort(pareto[,i]))
  rank<-c(rank,rank1)
  u.val<-c(u.val,u.values1)
}
df<-data.frame(rank=rank, u.values=u.val, 
               cell_expression=c(rep("low",
                                     length((unique((sort(pareto[,1])))))),
                                 rep("med", length((unique((sort(pareto[,2])))))),
                                 rep("high", length((unique((sort(pareto[,3]))))))
               ))
g12<-ggplot(df, aes(x=log1p(rank), y=log1p(u.values),
                    color=cell_expression))+
  geom_point(shape=1)+
  geom_smooth(method="lm", fullrange=T, se=FALSE)+ 
  scale_color_manual(values=c('#56B4E9','#E69F00', '#999999'))+
  theme_minimal(base_line_size = .75)+ 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

#CSEQ53
counts<-counts_3cel$res_cs53
rank<-c()
u.val<-c()
for (i in 1:ncol(counts)) {
  rank1<-rank(unique((sort(counts[,i], decreasing = T))))
  u.values1<-unique(sort(counts[,i]))
  rank<-c(rank,rank1)
  u.val<-c(u.val,u.values1)
}
df<-data.frame(rank=rank, u.values=u.val, 
               cell_expression=c(rep("low",
                                     length((unique((sort(counts[,1])))))),
                                 rep("med", length((unique((sort(counts[,2])))))),
                                 rep("high", length((unique((sort(counts[,3]))))))
               ))
g13<-ggplot(df, aes(x=log1p(rank), y=log1p(u.values),
                    color=cell_expression))+
  geom_point(shape=1)+
  geom_smooth(method="lm", fullrange=T, se=FALSE)+ 
  scale_color_manual(values=c('#56B4E9','#E69F00', '#999999'))+
  theme_minimal(base_line_size = .75)+ 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold",hjust = 0.5))+
  ggtitle("CELSeq53")

pareto<-par_3cel$res_cs53
rank<-c()
u.val<-c()
for (i in 1:ncol(pareto)) {
  rank1<-rank(unique((sort(pareto[,i], decreasing = T))))
  u.values1<-unique(sort(pareto[,i]))
  rank<-c(rank,rank1)
  u.val<-c(u.val,u.values1)
}
df<-data.frame(rank=rank, u.values=u.val, 
               cell_expression=c(rep("low",
                                     length((unique((sort(pareto[,1])))))),
                                 rep("med", length((unique((sort(pareto[,2])))))),
                                 rep("high", length((unique((sort(pareto[,3]))))))
               ))
g14<-ggplot(df, aes(x=log1p(rank), y=log1p(u.values),
                    color=cell_expression))+
  geom_point(shape=1)+
  geom_smooth(method="lm", fullrange=T, se=FALSE)+ 
  scale_color_manual(values=c('#56B4E9','#E69F00', '#999999'))+
  theme_minimal(base_line_size = .75)+ 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

#10x5
counts<-counts_3cel$res_10x5
rank<-c()
u.val<-c()
for (i in 1:ncol(counts)) {
  rank1<-rank(unique((sort(counts[,i], decreasing = T))))
  u.values1<-unique(sort(counts[,i]))
  rank<-c(rank,rank1)
  u.val<-c(u.val,u.values1)
}
df<-data.frame(rank=rank, u.values=u.val, 
               cell_expression=c(rep("low",
                                     length((unique((sort(counts[,1])))))),
                                 rep("med", length((unique((sort(counts[,2])))))),
                                 rep("high", length((unique((sort(counts[,3]))))))
               ))
g3<-ggplot(df, aes(x=log1p(rank), y=log1p(u.values),
                   color=cell_expression))+
  geom_point(shape=1)+
  geom_smooth(method="lm", fullrange=T, se=FALSE)+ 
  scale_color_manual(values=c('#56B4E9','#E69F00', '#999999'))+
  theme_minimal(base_line_size = .75)+ 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold",hjust = 0.5))+
  ggtitle("10x5")

pareto<-par_3cel$res_10x5
rank<-c()
u.val<-c()
for (i in 1:ncol(pareto)) {
  rank1<-rank(unique((sort(pareto[,i], decreasing = T))))
  u.values1<-unique(sort(pareto[,i]))
  rank<-c(rank,rank1)
  u.val<-c(u.val,u.values1)
}
df<-data.frame(rank=rank, u.values=u.val, 
               cell_expr=c(rep("Low Depth",
                               length((unique((sort(pareto[,1])))))),
                           rep("Med Depth", length((unique((sort(pareto[,2])))))),
                           rep("High Depth", length((unique((sort(pareto[,3]))))))
               ))
g4<-ggplot(df, aes(x=log1p(rank), y=log1p(u.values),
                   color=cell_expr))+
  geom_point(shape=1)+
  geom_smooth(method="lm", fullrange=T, se=FALSE)+ 
  scale_color_manual(values=c('#56B4E9','#E69F00', '#999999'))+
  theme_minimal(base_line_size = .75)+ 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
# ------------------- refactoring end ----------
library(grid)
library(ggplot2)
library(lattice)
library(gridExtra)
t1<-textGrob("Raw")
t2<-textGrob("PsiNorm")
r1 <- textGrob("CELSeq")
r2 <- textGrob("10x")
r3 <- textGrob("DropSeq")
r4 <- textGrob("CSeq1")
r5 <- textGrob("CSeq2")
r6 <- textGrob("CSeq3")
r7<-textGrob("10x5")
b1 <- textGrob("")

legend<-get_legend(ggplot(df, aes(x=log1p(rank), y=log1p(u.values),
                                  color=cell_expr))+
                     geom_point()+ 
                     scale_color_manual(values=c('#56B4E9',
                                                 '#E69F00', '#999999'))+
                     theme_minimal(base_line_size = .75)+ 
                     theme(legend.position = "right", 
                           legend.title = element_blank()))

ss<-grid.arrange(right=legend,t1,g1,g5,g7,
                 t2,g2,g6,g8, nrow=2, ncol=4,
                 bottom="Log(rank)", left="Log(expression)")
ss2<-grid.arrange(right=legend,t1,g9,g11,g13,g3,
                  t2,g10,g12,g14,g4, nrow=2, ncol=5,
                  bottom="Log(rank)", left="Log(expression)")
grid.arrange(ss, graf1) # Figure 2 for main text
grid.arrange(ss2, graf2) # for Supplementary

