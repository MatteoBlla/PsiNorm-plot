library(ggplot2)
library(magrittr)
library(CellBench)
library(SingleCellExperiment)
library(gridExtra)
library(VGAM)
library(cowplot)

#loading the datasets after the normalizations
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
  assay(datasets$res_10x, norm[i])<-assay(res1$result[[i+8]], norm [i])
  assay(datasets$res_dseq, norm[i])<-assay(res1$result[[i+16]], norm [i])
  assay(datasets$res_cs51, norm[i])<-assay(res1$result[[i+24]], norm [i])
  assay(datasets$res_cs52, norm[i])<-assay(res1$result[[i+32]], norm [i])
  assay(datasets$res_cs53, norm[i])<-assay(res1$result[[i+40]], norm [i])
}
remove(res1)
load("/directory/10x_5cl_after_normalization.rdata")
datasets$res_10x5<-res1$result[[1]]
for (i in 2:8) {
  assay(datasets$res_10x5, norm[i])<-assay(res1$result[[i]], norm [i])
}
remove(res1)

# counts matrix
datas<-lapply(datasets, function(x) counts(x)) 

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

rank_plot<-function(list){
  gg<-list()
  for (j in 1:length(list)){
    title<-c("CelSeq", "10x", "DropSeq", "CelSeq1",
             "CelSeq2","CelSeq3", "10x5")
    c<- list[[j]]
    u.val<-c()
    rank<-c()
    for (i in 1:ncol(c)){
      u.valu<-as.numeric(table(c[,i])) #how many time appears each value 
      ranku<-rank(unique((sort(c[,i], decreasing = F)))) #ranks
      u.val<-c(u.val,u.valu)
      rank<-c(rank, ranku)
    }
    df<-data.frame(rank=rank, val=u.val, 
                   cell_expression=c(rep("low",
                                         length(unique(sort(c[,1])))),
                                     rep("med", length(unique(sort(c[,2])))),
                                     rep("high", length(unique(sort(c[,3]))))
                   ))
    lm1<-lm(log(df$val[df$cell_expression=="low"])~log(df$rank[df$cell_expression=="low"])) #to compute the R^2
    lm2<-lm(log(df$val[df$cell_expression=="med"])~log(df$rank[df$cell_expression=="med"]))
    lm3<-lm(log(df$val[df$cell_expression=="high"])~log(df$rank[df$cell_expression=="high"]))
    gg[[j]]<-ggplot(df, aes(x=log(rank), y=log(val), 
                            color=cell_expression))+
      geom_point(shape=1)+
      ggtitle(title[j], 
              subtitle = paste0("slope: ",
                                round(summary(lm1)$coefficients[2], 
                                      digits = 3), 
                                " R^2: ", round(summary(lm1)$r.squared,
                                                digits = 3), " (L)\r\n",
                                "slope: ", round(summary(lm2)$coefficients[2], 
                                                 digits = 3), 
                                " R^2: ", round(summary(lm2)$r.squared,
                                                digits = 3), " (M)\r\n",
                                "slope: ", round(summary(lm3)$coefficients[2], 
                                                 digits = 3), 
                                " R^2: ", round(summary(lm3)$r.squared,
                                                digits = 3), " (H)"))+
      geom_smooth(method = lm, aes(fill=cell_expression), se=F)+
      theme_minimal(base_line_size = .75)+ 
      theme(legend.position = "none",
            axis.title.x = element_blank(), 
            axis.title.y = element_blank(),
            plot.title = element_text(size=20, face="bold",hjust = 0.5))+ 
      scale_color_manual(values=c('#56B4E9',
                                  '#E69F00', '#999999'))
    
  }
  gg
}

gg<-rank_plot(counts_3cel)

df<-data.frame(rank=gg[[7]]$data$rank, val=gg[[7]]$data$val, #to make the legend
               cell_expression=c(rep(" Low Depth",  
                                     length(unique(sort(counts_3cel$res_10x5[,1])))),
                                 rep(" Med Depth", length(unique(sort(counts_3cel$res_10x5[,2])))),
                                 rep("High Depth", length(unique(sort(counts_3cel$res_10x5[,3]))))
               ))
legend<-get_legend(ggplot(df, aes(x=log(rank), y=log(val),
                                  color=cell_expression))+
                     geom_point()+theme_minimal()+ theme(legend.title = element_blank())+
                     scale_color_manual(values=c('#E69F00',
                                                 '#999999', '#56B4E9')))
# FIRST PART OF THE FIGURE 1
slo1<-grid.arrange(#gg[[1]], gg[[2]], gg[[3]], right=legend, ncol=3, #comment one line at time
  gg[[4]] ,gg[[5]], gg[[6]], gg[[7]], right=legend, ncol=4, #for supplementary
  bottom="Log(ranks)",
  left="Log(absolute frequencies)",
  nrow=1)

# GOODNESS OF FIT
pareto.MLE <- function(X) #estimating Pareto factors
{
  n <- length(X[X>0])
  m <- min(X[X>0])
  a <- n/sum(log(X[X>0])-log(m))
  return(c(m,a)) 
}

# generating a Pareto from a Uniform
x.random <- function(N, x.min, alpha){
  x.min*(1-runif(N))^(-1/alpha)
}

# compute the differences between teorical quantiles and empirical ones
scarto.quantili<-function(data){
  N=nrow(data) 
  m=apply(data, 2, pareto.MLE)[1,] 
  alpha<-apply(data+1, 2, pareto.MLE)[2,] 
  alpha0<-apply(data, 2, pareto.MLE)[2,]
  s<-c()
  for(i in 1:ncol(data)){
    uno<-sort(data[,i], decreasing = T) 
    fr<-uno[uno>0]  
    p <- fr/sum(fr) 
    ll <- function(s)     
      sum(fr*(s*log(1:length(p))+log(sum(1/(1:length(p))^s)))) 
    fit <- mle(ll,start=list(s=1))
    s.ll <- coef(fit)
    s <- c(s,s.ll)
  }
  c<-data*0
  c0<-data*0
  c.zipf<-data*0
  
  q.pareto<-c() 
  q.pareto0<-c()
  q.zipf<-c()
  
  for (i in 1:ncol(data)){
    c[,i]<-round(x.random(N, m[i], alpha[i])) # pareto+1 simulation 
    c0[,i]<-round(x.random(N, m[i], alpha0[i])) # pareto0 simulation
    c.zipf[,i]<-rzipf(n=nrow(data), N=max(data[,i]), shape = s[i]) #simulation of zipf
    
    q.pareto[i]<-log((quantile(c[,i]-1, 0.75)+1)/(quantile(data[,i],0.75)+1))
    q.pareto0[i]<-log((quantile(c0[,i]-1, 0.75)+1)/(quantile(data[,i],0.75)+1))
    q.zipf[i]<-log((quantile(c.zipf[,i], 0.75)+1)/(quantile(data[,i],0.75)+1))
  }
  
  q.max<-cbind(q.pareto, q.pareto0, q.zipf)
  
  diff.quantili<-list(q.max)
  diff.quantili
}

scarti<-lapply(datas, scarto.quantili)

scarti.75<-list(scarti$res_cs2[[1]],
                scarti$res_10x[[1]],
                scarti$res_dseq[[1]],
                scarti$res_cs51[[1]],
                scarti$res_cs52[[1]],
                scarti$res_cs53[[1]],
                scarti$res_10x5[[1]])
names(scarti.75)<-names(scarti)

library(ggplot2)
library(ggridges)
#CELSEQ
df<-data.frame(scarti=as.vector(scarti.75$res_cs2),
               metodo=c(rep("Pareto+1", nrow(scarti.75$res_cs2)),
                        rep("Pareto 0", nrow(scarti.75$res_cs2)),
                        rep("Zipf", nrow(scarti.75$res_cs2))),
               ciao=rep("ciao", nrow(scarti.75$res_cs2)*3))
gcelseq<-ggplot(df, aes(x=scarti, y=metodo))+
  geom_density_ridges(alpha=.7, aes(fill=metodo))+
  theme_minimal(base_line_size = .75)+ 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())
#10X
df<-data.frame(scarti=as.vector(scarti.75$res_10x),
               metodo=c(rep("Pareto+1", nrow(scarti.75$res_10x)),
                        rep("Pareto 0", nrow(scarti.75$res_10x)),
                        rep("Zipf", nrow(scarti.75$res_10x))),
               ciao=rep("ciao", nrow(scarti.75$res_10x)*3))
g10x<-ggplot(df, aes(x=scarti, y=metodo))+
  geom_density_ridges(alpha=.7, aes(fill=metodo))+
  theme_minimal(base_line_size = .75)+ 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

#DROPSEQ 
df<-data.frame(scarti=as.vector(scarti.75$res_dseq),
               metodo=c(rep("Pareto+1", nrow(scarti.75$res_dseq)),
                        rep("Pareto 0", nrow(scarti.75$res_dseq)),
                        rep("Zipf", nrow(scarti.75$res_dseq))),
               ciao=rep("ciao", nrow(scarti.75$res_dseq)*3))
gdrop<-ggplot(df, aes(x=scarti, y=metodo))+
  geom_density_ridges(alpha=.6, aes(fill=metodo))+
  theme_minimal(base_line_size = .75)+ 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

#CSEQ1
df<-data.frame(scarti=as.vector(scarti.75$res_cs51),
               metodo=c(rep("Pareto+1", nrow(scarti.75$res_cs51)),
                        rep("Pareto 0", nrow(scarti.75$res_cs51)),
                        rep("Zipf", nrow(scarti.75$res_cs51))),
               ciao=rep("ciao", nrow(scarti.75$res_cs51)*3))
gcs1<-ggplot(df, aes(x=scarti, y=metodo))+
  geom_density_ridges(alpha=.6, aes(fill=metodo))+
  theme_minimal(base_line_size = .75)+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())
#CSEQ2
df<-data.frame(scarti=as.vector(scarti.75$res_cs52),
               metodo=c(rep("Pareto+1", nrow(scarti.75$res_cs52)),
                        rep("Pareto 0", nrow(scarti.75$res_cs52)),
                        rep("Zipf", nrow(scarti.75$res_cs52))),
               ciao=rep("ciao", nrow(scarti.75$res_cs52)*3))
gcs2<-ggplot(df, aes(x=scarti, y=metodo))+
  geom_density_ridges(alpha=.6, aes(fill=metodo))+
  theme_minimal(base_line_size = .75)+ 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())
#CSEQ3
df<-data.frame(scarti=as.vector(scarti.75$res_cs53),
               metodo=c(rep("Pareto+1", nrow(scarti.75$res_cs53)),
                        rep("Pareto 0", nrow(scarti.75$res_cs53)),
                        rep("Zipf", nrow(scarti.75$res_cs53))),
               ciao=rep("ciao", nrow(scarti.75$res_cs53)*3))
gcs3<-ggplot(df, aes(x=scarti, y=metodo))+
  geom_density_ridges(alpha=.6, aes(fill=metodo))+
  theme_minimal(base_line_size = .75)+ 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())
#10X5
df<-data.frame(scarti=as.vector(scarti.75$res_10x5),
               metodo=c(rep("Pareto+1", nrow(scarti.75$res_10x5)),
                        rep("Pareto 0", nrow(scarti.75$res_10x5)),
                        rep("Zipf", nrow(scarti.75$res_10x5))),
               ciao=rep("ciao", nrow(scarti.75$res_10x5)*3))
g10x5<-ggplot(df, aes(x=scarti, y=metodo))+
  geom_density_ridges(alpha=.6, aes(fill=metodo))+
  theme_minimal(base_line_size = .75)+ 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

legend<-get_legend(ggplot(df, aes(x=scarti, y=metodo))+
                     geom_density_ridges(alpha=.6, aes(fill=metodo),
                                         stat = "binline", bins = 50)+
                     theme(legend.title = element_blank()))  

g1<-grid.arrange(#gcelseq, g10x, gdrop, ncol=3, #comment each row at time
  gcs1, gcs2, gcs3, g10x5, ncol=4,  #for supplementary
  bottom="Log ratio simulated versus empirical third quartile",
  left="Absolute Frequencies",
  right=legend)    
# FIGURE 1
plot_grid(slo1,g1, ncol=1, labels = c("A", "B")) 
