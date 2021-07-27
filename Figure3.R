### 25k cells
datas<-list(
  scran=read.csv(file = "/directory/scran25.txt", sep = ""),
  pareto=read.csv(file = "/directory/pareto25.txt", sep = ""),
  deseq=read.csv(file = "/directory/deseq25.txt", sep = ""),
  logcpm=read.csv(file = "/directory/logcpm25.txt", sep = ""),
  clr=read.csv(file = "/directory/clr25.txt", sep = ""),
  linnorm=read.csv(file = "/directory/linnorm25.txt", sep = ""),
  tmm=read.csv(file = "/directory/tmm25.txt", sep = ""),
  sct=read.csv(file = "/directory/NEWsct25.txt", sep = "")
)

datas<-lapply(datas, function(x) x[which(x$User=="matteo"),])
lapply(datas, function(x) table(x$PID)) 
datas$logcpm<-datas$logcpm[datas$logcpm$PID==30031,]

library(stringr)
numextract <- function(string){ 
  str_extract(string, "\\-*\\d+\\.*\\d*")
}


datas<-lapply(datas, function(x) numextract(x$RSS))
datas<-lapply(datas, function(x) as.numeric(x))
pippo<-lapply(datas, function(x) x=c(x[x>100]/1000, x[x<100]))
pippo<-lapply(pippo, function(x) sort(x, decreasing = F))

df<-data.frame(RAM=c(pippo$scran, pippo$pareto, pippo$deseq, pippo$logcpm,
                     pippo$clr, pippo$linnorm, pippo$tmm, pippo$sct),
               norm=c(rep("Scran", length(pippo$scran)),
                      rep("PsiNorm", length(pippo$pareto)),
                      rep("DESeq", length(pippo$deseq)),
                      rep("logCPM", length(pippo$logcpm)),
                      rep("CLR", length(pippo$clr)),
                      rep("Linnorm", length(pippo$linnorm)),
                      rep("TMM", length(pippo$tmm)),
                      rep("sctransform", length(pippo$sct))),
               time=c(seq(1:length(pippo$scran)),
                      seq(1:length(pippo$pareto)),
                      seq(1:length(pippo$deseq)),
                      seq(1:length(pippo$logcpm)),
                      seq(1:length(pippo$clr)),
                      seq(1:length(pippo$linnorm)),
                      seq(1:length(pippo$tmm)),
                      seq(1:length(pippo$sct))
               ))
df$norm<-as.factor(df$norm)
ram25<-rep(NA,8)
tim25<-rep(NA,8)
names(ram25)<-levels(df$norm)
names(tim25)<-levels(df$norm)
for (i in 1:8) {
  ram25[i]<-max(df$RAM[df$norm==levels(df$norm)[i]])
  tim25[i]<-max(df$time[df$norm==levels(df$norm)[i]])
}

scran10=read.csv(file = "/directory/scran25_10c.txt", sep = "") #scran with 10 core
strptime(scran10$X.Time, format = "%H:%M:%S")[577]-strptime(scran10$X.Time, format = "%H:%M:%S")[1]

tim25<-c(tim25, 2.9*60)
names(tim25)<-c(names(tim25)[1:8], "Scran10")
scran10$RSS<-as.numeric(numextract(scran10$RSS))
scran10$RSS[scran10$RSS>100]<-scran10$RSS[scran10$RSS>100]/1000
ram25<-c(ram25, max(scran10$RSS))
names(ram25)<-c(names(ram25)[1:8], "Scran10")

### 50k
datas<-list(
  scran=read.csv(file = "/directory/scran50.txt", sep = ""),
  pareto=read.csv(file = "/directory/pareto50.txt", sep = ""),
  deseq=read.csv(file = "/directory/deseq50.txt", sep = ""),
  logcpm=read.csv(file = "/directory/logcpm50.txt", sep = ""),
  clr=read.csv(file = "/directory/clr50.txt", sep = ""),
  linnorm=read.csv(file = "/directory/linnorm50.txt", sep = ""),
  tmm=read.csv(file = "/directory/tmm50.txt", sep = ""),
  sct=read.csv(file = "/directory/NEWsct50.txt", sep = "")
)

datas<-lapply(datas, function(x) x[which(x$User=="matteo"),])
lapply(datas, function(x) table(x$PID)) 

datas<-lapply(datas, function(x) numextract(x$RSS))
datas<-lapply(datas, function(x) as.numeric(x))
pippo<-lapply(datas, function(x) x=c(x[x>100]/1000, x[x<100]))
pippo<-lapply(pippo, function(x) sort(x, decreasing = F))

df<-data.frame(RAM=c(pippo$scran, pippo$pareto, pippo$deseq, pippo$logcpm,
                     pippo$clr, pippo$linnorm, pippo$tmm, pippo$sct),
               norm=c(rep("Scran", length(pippo$scran)),
                      rep("PsiNorm", length(pippo$pareto)),
                      rep("DESeq", length(pippo$deseq)),
                      rep("logCPM", length(pippo$logcpm)),
                      rep("CLR", length(pippo$clr)),
                      rep("Linnorm", length(pippo$linnorm)),
                      rep("TMM", length(pippo$tmm)),
                      rep("sctransform", length(pippo$sct))),
               time=c(seq(1:length(pippo$scran)),
                      seq(1:length(pippo$pareto)),
                      seq(1:length(pippo$deseq)),
                      seq(1:length(pippo$logcpm)),
                      seq(1:length(pippo$clr)),
                      seq(1:length(pippo$linnorm)),
                      seq(1:length(pippo$tmm)),
                      seq(1:length(pippo$sct))
               ))
df$norm<-as.factor(df$norm)

ram50<-rep(NA,8)
tim50<-rep(NA,8)
names(ram50)<-levels(df$norm)
names(tim50)<-levels(df$norm)
for (i in 1:8) {
  ram50[i]<-max(df$RAM[df$norm==levels(df$norm)[i]])
  tim50[i]<-max(df$time[df$norm==levels(df$norm)[i]])
}

scran10=read.csv(file = "/directory/scran50_10c.txt", sep = "") #scran 10 core
strptime(scran10$X.Time, format = "%H:%M:%S")[length(scran10$X.Time)]-strptime(scran10$X.Time, format = "%H:%M:%S")[1]

tim50<-c(tim50, 5.133333*60)
names(tim50)<-c(names(tim50)[1:8], "Scran10")
scran10$RSS<-as.numeric(numextract(scran10$RSS))
scran10$RSS[scran10$RSS>100]<-scran10$RSS[scran10$RSS>100]/1000
ram50<-c(ram50, max(scran10$RSS))
names(ram50)<-c(names(ram50)[1:8], "Scran10")


#### 75 000 cells 
datas<-list(
  scran=read.csv(file = "/directory/scran75.txt", sep = ""),
  pareto=read.csv(file = "/directory/pareto75.txt", sep = ""),
  deseq=read.csv(file = "/directory/deseq75.txt", sep = ""),
  logcpm=read.csv(file = "/directory/logcpm75.txt", sep = ""),
  clr=read.csv(file = "/directory/clr75.txt", sep = ""),
  linnorm=read.csv(file = "/directory/linnorm75.txt", sep = ""),
  tmm=read.csv(file = "/directory/tmm75.txt", sep = ""),
  sct=read.csv(file = "/directory/NEWsct75.txt", sep = "")
)

datas<-lapply(datas, function(x) x[which(x$User=="matteo"),])
lapply(datas, function(x) table(x$PID)) 
datas$deseq<-datas$deseq[datas$deseq$PID==29777,]

datas<-lapply(datas, function(x) numextract(x$RSS))
datas<-lapply(datas, function(x) as.numeric(x))
pippo<-lapply(datas, function(x) x=c(x[x>100]/1000, x[x<100]))
pippo<-lapply(pippo, function(x) sort(x, decreasing = F))

df<-data.frame(RAM=c(pippo$scran, pippo$pareto, pippo$deseq, pippo$logcpm,
                     pippo$clr, pippo$linnorm, pippo$tmm, pippo$sct),
               norm=c(rep("Scran", length(pippo$scran)),
                      rep("PsiNorm", length(pippo$pareto)),
                      rep("DESeq", length(pippo$deseq)),
                      rep("logCPM", length(pippo$logcpm)),
                      rep("CLR", length(pippo$clr)),
                      rep("Linnorm", length(pippo$linnorm)),
                      rep("TMM", length(pippo$tmm)),
                      rep("sctransform", length(pippo$sct))),
               time=c(seq(1:length(pippo$scran)),
                      seq(1:length(pippo$pareto)),
                      seq(1:length(pippo$deseq)),
                      seq(1:length(pippo$logcpm)),
                      seq(1:length(pippo$clr)),
                      seq(1:length(pippo$linnorm)),
                      seq(1:length(pippo$tmm)),
                      seq(1:length(pippo$sct))
               ))
df$norm<-as.factor(df$norm)

ram75<-rep(NA,8)
tim75<-rep(NA,8)
names(ram75)<-levels(df$norm)
names(tim75)<-levels(df$norm)
for (i in 1:8) {
  ram75[i]<-max(df$RAM[df$norm==levels(df$norm)[i]])
  tim75[i]<-max(df$time[df$norm==levels(df$norm)[i]])
}

scran10=read.csv(file = "/directory/scran75_10c.txt", sep = "") #scran with 10 core
strptime(scran10$X.Time, format = "%H:%M:%S")[length(scran10$X.Time)]-strptime(scran10$X.Time, format = "%H:%M:%S")[1]

tim75<-c(tim75, 7.4*60)
names(tim75)<-c(names(tim75)[1:8], "Scran10")
scran10$RSS<-as.numeric(numextract(scran10$RSS))
scran10$RSS[scran10$RSS>100]<-scran10$RSS[scran10$RSS>100]/1000
ram75<-c(ram75, max(scran10$RSS))
names(ram75)<-c(names(ram75)[1:8], "Scran10")

## 100
datas<-list(
  scran=read.csv(file = "/directory/scran100.txt", sep = ""),
  pareto=read.csv(file = "/directory/pareto100.txt", sep = ""),
  deseq=read.csv(file = "/directory/deseq100.txt", sep = ""),
  logcpm=read.csv(file = "/directory/logcpm100.txt", sep = ""),
  clr=read.csv(file = "/directory/clr100.txt", sep = ""),
  linnorm=read.csv(file = "/directory/linnorm100.txt", sep = ""),
  tmm=read.csv(file = "/directory/tmm100.txt", sep = ""),
  sct=read.csv(file = "/directory/NEWsct100.txt", sep = "")
)

datas<-lapply(datas, function(x) x[which(x$User=="matteo"),])
lapply(datas, function(x) table(x$PID)) 
datas$linnorm<-datas$linnorm[datas$linnorm$PID==31088,]

datas<-lapply(datas, function(x) numextract(x$RSS))
datas<-lapply(datas, function(x) as.numeric(x))
pippo<-lapply(datas, function(x) x=c(x[x>100]/1000, x[x<100]))
pippo<-lapply(pippo, function(x) sort(x, decreasing = F))

df<-data.frame(RAM=c(pippo$scran, pippo$pareto, pippo$deseq, pippo$logcpm,
                     pippo$clr, pippo$linnorm, pippo$tmm, pippo$sct),
               norm=c(rep("Scran", length(pippo$scran)),
                      rep("PsiNorm", length(pippo$pareto)),
                      rep("DESeq", length(pippo$deseq)),
                      rep("logCPM", length(pippo$logcpm)),
                      rep("CLR", length(pippo$clr)),
                      rep("Linnorm", length(pippo$linnorm)),
                      rep("TMM", length(pippo$tmm)),
                      rep("sctransform", length(pippo$sct))),
               time=c(seq(1:length(pippo$scran)),
                      seq(1:length(pippo$pareto)),
                      seq(1:length(pippo$deseq)),
                      seq(1:length(pippo$logcpm)),
                      seq(1:length(pippo$clr)),
                      seq(1:length(pippo$linnorm)),
                      seq(1:length(pippo$tmm)),
                      seq(1:length(pippo$sct))
               ))
df$norm<-as.factor(df$norm)

ram100<-rep(NA,8)
tim100<-rep(NA,8)
names(ram100)<-levels(df$norm)
names(tim100)<-levels(df$norm)
for (i in 1:8) {
  ram100[i]<-max(df$RAM[df$norm==levels(df$norm)[i]])
  tim100[i]<-max(df$time[df$norm==levels(df$norm)[i]])
}

scran10=read.csv(file = "/directory/scran100_10c.txt", sep = "")
strptime(scran10$X.Time, format = "%H:%M:%S")[length(scran10$X.Time)]-strptime(scran10$X.Time, format = "%H:%M:%S")[1]

tim100<-c(tim100, 10.5333333*60)
names(tim100)<-c(names(tim100)[1:8], "Scran10")
scran10$RSS<-as.numeric(numextract(scran10$RSS))
scran10$RSS[scran10$RSS>100]<-scran10$RSS[scran10$RSS>100]/1000
ram100<-c(ram100, max(scran10$RSS))
names(ram100)<-c(names(ram100)[1:8], "Scran10")


df<-data.frame(RAM=c(ram25, ram50, ram75, ram100),
               time=c(tim25, tim50, tim75, tim100),
               norm=rep(names(ram25), 4),
               ncell=c(rep(" 25.000",9),
                       rep(" 50.000",9),
                       rep(" 75.000",9),
                       rep("100.000",9)))
df$norm[df$norm=="Scran10"]<-"Scran (10 core)"

library("RColorBrewer")
library(ggplot2)
gg1<-ggplot(df, aes(x=2*time, y=RAM))+
  geom_point(aes(color=norm, size=ncell))+
  scale_colour_manual(values=c("#660000","#FF6600",
                               "#66FF33","#009900",
                               "#3399FF","#000099",
                               "#9900FF","#00CCFF", "666666"))+
  geom_line(aes(color=norm))+
  labs(x="Time(sec)", y="RAM(Gb)")+ 
  scale_x_continuous(breaks = 240*c(0,1,2,3,4,5,6,7,8,9))+
  theme_minimal(base_line_size = .75)+
  theme(legend.position = "none",
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size=13))

a<-which(df$norm=="PsiNorm")
b<-which(df$norm=="logCPM")
c<-which(df$norm=="DESeq")
d<-which(df$norm=="Linnorm")
e<-which(df$norm=="TMM")
f<-which(df$norm=="CLR")
g<-which(df$norm=="sctransform")
df2<-df[c(a,b,c,d,e,f,g),]
ggplot(df2, aes(x=2*time, y=RAM))+
  geom_point(aes(color=norm, size=ncell))+
  scale_colour_manual(values=c("#660000","#FF6600",
                               "#66FF33","#009900",
                               "#3399FF","#00CCFF","666666"))+
  geom_line(aes(color=norm))+
  scale_x_continuous(breaks = 240*c(0,1,2,3,4,5,6,7,8,9))+
  theme_minimal(base_line_size = .75)+ 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        legend.position = "none")

#time by number of cells
df<-data.frame(time=c(tim25, tim50, tim75, tim100),
               norm=rep(names(tim25), times=4),
               ncell=rep(c(25000,50000,75000,100000), each=9))

gtime<-ggplot(df, aes(ncell, time))+
  geom_point(aes(color=norm, size=ncell))+
  theme_minimal()+
  geom_line(aes(color=norm))+
  theme(legend.position = "none")+
  scale_colour_manual(values=c("#660000","#FF6600",
                               "#66FF33","#009900",
                               "#3399FF","#000099",
                               "#9900FF","#00CCFF", "666666"))+
  scale_x_continuous(breaks = c(25000,50000,75000,100000))+
  labs(x="Number of cells", y="Time (sec)")+
  theme(legend.position = "none",
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size=13))

library(cowplot)
library(gridExtra)
df<-data.frame(time=c(tim25, tim50, tim75, tim100),
               norm=rep(names(tim25), times=4),
               ncell=rep(c(" 25.000", " 50.000", " 75.000", "100.000"), each=9))
legend<-get_legend(ggplot(df, aes(ncell, time))+
                     geom_point(aes(color=norm, size=as.factor(ncell)))+
                     theme_minimal()+
                     theme(legend.title = element_blank(),
                           legend.text = element_text(size = 11),
                           legend.position = "bottom")+
                     geom_line(aes(color=norm))+
                     scale_colour_manual(values=c("#660000","#FF6600",
                                                  "#66FF33","#009900",
                                                  "#3399FF","#000099",
                                                  "#9900FF","#00CCFF", "666666")))
pg<-plot_grid(gg1, gtime, ncol=2, align = "hv", labels = c("A","B"))
grid.arrange(pg, bottom=legend)
