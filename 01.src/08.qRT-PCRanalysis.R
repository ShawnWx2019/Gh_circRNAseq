#################################
#   Prj: Rscript
#   Assignment: qPCR analysis
#   Date: Nov 30, 2020
#   Author: Shawn Wang
#################################
##=====Notice==========
## 01.适用于一种重复数据（生物学重复或者技术重复）
## 02.参考y叔远古时期的一篇博文，这里对于差异显著性计算用wilcoxon Rank Sum test
library(getopt)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(tidyverse)
library(cowplot)
library(ggpol)
library(grid)
library(rlang)
filter(data,Primer == "circ_2010")
ddeltaCt = function(data,reference, gene, material, Repl){
  ## ΔCt1: gene Ct value - reference Ct value.
  ## ΔCt2: average of ΔCt1 of control group.
  RG = list()
  TG = list()
  r = list()
  t = list()
  i = 1
  deltaCt1 = list()
  for (i in 1:length(material)) {
    RG[[i]] = filter(data,Primer == reference, Material == material[i])
    r[[i]] = mean(as.numeric(RG[[i]][,3]))
    TG[[i]] = filter(data,Primer == gene, Material == material[i])
    t[[i]] = mean(as.numeric(TG[[i]][,3]))
    ## ΔCt1 
    deltaCt1[[i]] = as.numeric(TG[[i]][,3]) - as.numeric(RG[[i]][,3])
  }
  ## average ct value of control group 
  ## gene
  r = r[[1]]
  meanr = mean(r)
  ## control
  t = t[[1]]
  meant = mean(t)
  deltaCt2 = meant - meanr
  ## delta delta ct
  deltaCt1 = unlist(deltaCt1)
  expmat = data.frame(Material = factor(rep(material,each = Repl),levels = material),
                      deltaCt1 = deltaCt1,
                      deltaCt2 = deltaCt2)
  expmat$NegddeltaCt = -(expmat$deltaCt1-expmat$deltaCt2)
  expmat$relativeExp = 2**(expmat$NegddeltaCt)
  n = 1:length(material)
  m = list()
  for (i in 1:(length(n)-1)) {
    m[[i]] = n[-c(1:i)]
  }
  m = unlist(m)
  l = rep(c(1:(length(n)-1)),times = c((length(n)-1):1))
  compl = list()
  for (i in 1:length(m)) {
    compl[[i]] = c(l[i],m[i]) 
  }
  p1 = ggplot(expmat,aes(x = Material,y = relativeExp))+
    geom_boxplot()+
    geom_jitter()+
    geom_signif(comparisons = compl,map_signif_level =T,test = t.test)+
    theme_bw()+
    theme(title = element_text(size = 16),
          legend.position = "none",
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14))+
    labs(x = "",y = expression(2-ΔΔct))
  
  p2 = ggplot(expmat)+
    geom_boxjitter(
      aes( x = Material,y = relativeExp, fill = Material),
      jitter.shape = 21, jitter.color = NA, jitter.size = 1, 
      jitter.params = list(height = 0, width = 0.04),
      outlier.color = NA, errorbar.draw = TRUE
    )+
    theme_bw()+
    theme(title = element_text(size = 16),
          legend.position = "none",
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14))+
    labs(x = "",y = expression(2-ΔΔct))
  a = list(expmat = expmat,
           plot1 = p1,
           plot2 = p2)
  return(a)
}
## import data
x = read.table("~/03.project/03.circRNA/CircRNAVersion3/02.data/03.workingdir/result3.txt",
               sep = "\t",
               header = T)
#View(x)
reference = "UBQ"
material = c("BTF","BTJF","JF")
Repl = 3
gene = "circ_0968"
y =ddeltaCt(data= x,
            reference = reference,
            gene = gene,
            material = material,
            Repl = Repl)

y$expmat
y$plot1
y$plot2
p1 = y$plot1
p2 = y$plot2
ggsave(p1,filename = "/Volumes/Samsung_T5/03.circRNA/CircRNAVersion3/03.process/08.qPCR/0051_O.pdf",width = 4,height = 4,dpi = 300)
ggsave(p2,filename = "/Volumes/Samsung_T5//03.circRNA/CircRNAVersion3/03.process/08.qPCR/0051_O.pdf",width = 4,height = 4,dpi = 300)


x1 = read.table("~/03.project/03.circRNA/CircRNAVersion3/03.process/08.qPCR/qPCR-1.txt",
                header = T,
                sep = "\t",
                stringsAsFactors = F)
head(x1)
#View(x1)
names(x1) = c("Material","CT","Primer")
x1 = x1[,c(1,3,2)]
reference = "UBQ"
material = c("Bt-F","BtJ-F","J-F")
Repl = 3
gene = "circ_0968"
y =ddeltaCt(data = x1,
            reference = reference,
            gene = gene,
            material = material,
            Repl = Repl)
y$expmat
y$plot1
y$plot2
p1 = y$plot1
p2 = y$plot2
ggsave(p1,filename = "~/03.project/03.circRNA/CircRNAVersion3/03.process/08.qPCR/0968_F.pdf",width = 4,height = 4,dpi = 300)
ggsave(p2,filename = "~/03.project/03.circRNA/CircRNAVersion3/03.process/08.qPCR/bj1639-2_F.pdf",width = 4,height = 4,dpi = 300)

x2 = read.table("~/03.project/03.circRNA/CircRNAVersion3/03.process/08.qPCR/qPCR-2.txt",
                header = T,
                sep = "\t",
                stringsAsFactors = F)
head(x2)
names(x2) = c("CT","Material","Primer")
x2 = x2[,c(2,3,1)]
reference = "UBQ"
material = c("Bt-L","BtJ-L","J-L")
Repl = 3
gene = "circ-1194-2"
y =ddeltaCt(data = x2,
            reference = reference,
            gene = gene,
            material = material,
            Repl = Repl)
y$expmat
y$plot1
y$plot2
p1 = y$plot1
p2 = y$plot2
ggsave(p1,filename = "~/03.project/03.circRNA/CircRNAVersion3/03.process/08.qPCR/1194-2_L.pdf",width = 4,height = 4,dpi = 300)
ggsave(p2,filename = "~/03.project/03.circRNA/CircRNAVersion3/03.process/08.qPCR/bj1194-2_F.pdf",width = 4,height = 4,dpi = 300)
unique(x2$Primer)
#corr of RNAseq and qRT-PCR
# RBtvsBtJ = log2(4326/1)
# RJvsBtJ = log2(1/1)
# RBtvsJ = log2(4326/1)
# qBtvsBtJ = log2((mean(y$expmat[,5][4:6])+1)/(mean(y$expmat[,5][1:3])+1))
# qJvsBtJ = log2((mean(y$expmat[,5][4:6])+1)/(mean(y$expmat[,5][7:9])+1))
# qBtvsJ = log2((mean(y$expmat[,5][1:3])+1)/(mean(y$expmat[,5][7:9])+1))
# Expvalue = data.frame(type = c("BtJvsBt","BtJvsJ","BtvsJ"),
#                       RNAseq = c(RBtvsBtJ,RJvsBtJ,RBtvsJ),
#                       qPCR = c(qBtvsBtJ,qJvsBtJ,qBtvsJ))
# lm_eqn = function(df){
#   m = lm(df[,2] ~ df[,3] , df);
#   eq <- substitute(italic(R)^2~"="~r2~","~italic(P)~"="~p,
#                    list(a = format(as.numeric(coef(m)[1]), digits = 2),
#                         b = format(as.numeric(coef(m)[2]), digits = 2),
#                         r2 = format(summary(m)$r.squared, digits = 3),
#                         p = format(summary(m)$coefficients[2,4], digits = 4)))
#   as.character(as.expression(eq));
# }
# lm_eqn(Expvalue)
# ggplot(data = Expvalue,mapping = aes(x = RNAseq, y = qPCR))+
#   geom_point(size = 1)+
#   geom_smooth(method="lm")+
#   geom_text(aes(x=-2,y=0.7,label=lm_eqn(Expvalue)),parse=T)+
#   theme_light()+
#   geom_hline(aes(yintercept=0), colour="black", linetype="dashed")+
#   geom_vline(aes(xintercept = 0), color = "black", linetype = "dashed")+
#   xlab(expression(paste("RNAseq  ",log[2],"(foldchange)")))+
#   ylab(expression(paste("mRNA  ",log[2],"(foldchange)")))

# perm.test <- function(d1, d2, nIter=10000) {
#   if (length(d2) == 1)
#     d2 <- rep(d2, length(d1))
#   m <- mean(d2-d1)
#   pooledData <- c(d1, d2)
#   n <- length(d1)
#   meanDiff <- numeric(nIter)
#   for (i in 1:nIter) {
#     idx <- sample(1:length(pooledData), n, replace=FALSE)
#     d1 <- pooledData[idx]
#     d2 <- pooledData[-idx]
#     meanDiff[i] <- mean(d2) - mean(d1)
#   }
#   p <- mean(abs(meanDiff) >= abs(m))
#   return(p)
# }