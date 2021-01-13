### 25k
datas<-list(
  scran=read.csv(file = "C:/Users/Matteo/Desktop/final_box/RAM25/RAM/scran25.txt", sep = ""),
  pareto=read.csv(file = "C:/Users/Matteo/Desktop/final_box/RAM25/RAM/pareto25.txt", sep = ""),
  deseq=read.csv(file = "C:/Users/Matteo/Desktop/final_box/RAM25/RAM/deseq25.txt", sep = ""),
  logcpm=read.csv(file = "C:/Users/Matteo/Desktop/final_box/RAM25/RAM/logcpm25.txt", sep = ""),
  clr=read.csv(file = "C:/Users/Matteo/Desktop/final_box/RAM25/RAM/clr25.txt", sep = ""),
  linnorm=read.csv(file = "C:/Users/Matteo/Desktop/final_box/RAM25/RAM/linnorm25.txt", sep = ""),
  tmm=read.csv(file = "C:/Users/Matteo/Desktop/final_box/RAM25/RAM/tmm25.txt", sep = ""),
  sct=read.csv(file = "C:/Users/Matteo/Desktop/final_box/RAM25/RAM/sct25.txt", sep = "")
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


### 50k
datas<-list(
  scran=read.csv(file = "C:/Users/Matteo/Desktop/final_box/RAM25/RAM/scran50.txt", sep = ""),
  pareto=read.csv(file = "C:/Users/Matteo/Desktop/final_box/RAM25/RAM/pareto50.txt", sep = ""),
  deseq=read.csv(file = "C:/Users/Matteo/Desktop/final_box/RAM25/RAM/deseq50.txt", sep = ""),
  logcpm=read.csv(file = "C:/Users/Matteo/Desktop/final_box/RAM25/RAM/logcpm50.txt", sep = ""),
  clr=read.csv(file = "C:/Users/Matteo/Desktop/final_box/RAM25/RAM/clr50.txt", sep = ""),
  linnorm=read.csv(file = "C:/Users/Matteo/Desktop/final_box/RAM25/RAM/linnorm50.txt", sep = ""),
  tmm=read.csv(file = "C:/Users/Matteo/Desktop/final_box/RAM25/RAM/tmm50.txt", sep = ""),
  sct=read.csv(file = "C:/Users/Matteo/Desktop/final_box/RAM25/RAM/sct50.txt", sep = "")
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


datas<-list(
  scran=read.csv(file = "C:/Users/Matteo/Desktop/final_box/RAM25/RAM/scran75.txt", sep = ""),
  pareto=read.csv(file = "C:/Users/Matteo/Desktop/final_box/RAM25/RAM/pareto75.txt", sep = ""),
  deseq=read.csv(file = "C:/Users/Matteo/Desktop/final_box/RAM25/RAM/deseq75.txt", sep = ""),
  logcpm=read.csv(file = "C:/Users/Matteo/Desktop/final_box/RAM25/RAM/logcpm75.txt", sep = ""),
  clr=read.csv(file = "C:/Users/Matteo/Desktop/final_box/RAM25/RAM/clr75.txt", sep = ""),
  linnorm=read.csv(file = "C:/Users/Matteo/Desktop/final_box/RAM25/RAM/linnorm75.txt", sep = ""),
  tmm=read.csv(file = "C:/Users/Matteo/Desktop/final_box/RAM25/RAM/tmm75.txt", sep = ""),
  sct=read.csv(file = "C:/Users/Matteo/Desktop/final_box/RAM25/RAM/sct75.txt", sep = "")
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

datas<-list(
  scran=read.csv(file = "C:/Users/Matteo/Desktop/final_box/RAM25/RAM/scran100.txt", sep = ""),
  pareto=read.csv(file = "C:/Users/Matteo/Desktop/final_box/RAM25/RAM/pareto100.txt", sep = ""),
  deseq=read.csv(file = "C:/Users/Matteo/Desktop/final_box/RAM25/RAM/deseq100.txt", sep = ""),
  logcpm=read.csv(file = "C:/Users/Matteo/Desktop/final_box/RAM25/RAM/logcpm100.txt", sep = ""),
  clr=read.csv(file = "C:/Users/Matteo/Desktop/final_box/RAM25/RAM/clr100.txt", sep = ""),
  linnorm=read.csv(file = "C:/Users/Matteo/Desktop/final_box/RAM25/RAM/linnorm100.txt", sep = ""),
  tmm=read.csv(file = "C:/Users/Matteo/Desktop/final_box/RAM25/RAM/tmm100.txt", sep = ""),
  sct=read.csv(file = "C:/Users/Matteo/Desktop/final_box/RAM25/RAM/sct100.txt", sep = "")
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

df<-data.frame(RAM=c(ram25,ram50,ram75,ram100),
               time=c(tim25,tim50,tim75,tim100),
               norm=rep(names(ram25), 8),
               ncell=c(rep(" 25mila",8),
                       rep(" 50mila",8),
                       rep(" 75mila",8),
                       rep("100mila",8)))

library("RColorBrewer")
library(ggplot2)
ggplot(df, aes(x=2*time, y=RAM))+
  geom_point(aes(color=norm, size=ncell))+
  scale_colour_manual(values=c("#660000","#FF6600",
                               "#66FF33","#009900",
                               "#3399FF","#000099",
                               "#9900FF","#00CCFF"))+
  geom_line(aes(color=norm))+
  labs(x="Time(sec)", y="RAM(Gb)")+ 
  scale_x_continuous(breaks = 240*c(0,1,2,3,4,5,6,7,8,9))+
  theme_minimal(base_line_size = .75)

a<-which(df$norm=="PsiNorm")
b<-which(df$norm=="logCPM")
c<-which(df$norm=="DESeq")
d<-which(df$norm=="Linnorm")
e<-which(df$norm=="TMM")
f<-which(df$norm=="CLR")
df2<-df[c(a,b,c,d,e,f),]
ggplot(df2, aes(x=2*time, y=RAM))+
  geom_point(aes(color=norm, size=ncell))+
  scale_colour_manual(values=c("#660000","#FF6600",
                               "#66FF33","#009900",
                               "#3399FF","#00CCFF"))+
  geom_line(aes(color=norm))+
  scale_x_continuous(breaks = 240*c(0,1,2,3,4,5,6,7,8,9))+
  theme_minimal(base_line_size = .75)+ 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        legend.position = "none")

