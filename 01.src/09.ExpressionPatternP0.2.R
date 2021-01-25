#################################################
#       Prj: Heterosis-circRNA
#       Assignment: Different expression analysis
#       Author: Shawn Wang
#       Date: Sep 04, 2020
#       Reference: https://github.com/Wendellab/CisTransRegulation/blob/master/inheritanceMode.r
#       Reference DOI: | https://doi.org/10.1038/s41467-019-13386-w 
#       Reference Author: Guanjing Hu, Ying Bao & Jonathan F. Wendel
#################################################
## notation: readcount value was update in Sep 01, 2020; the raw 
setwd("~/03.project/03.circRNA/CircRNAVersion3/03.process/02.ExpresTest")
## prepare
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(svglite)
library(ggpubr)
options(stringsAsFactors = F)

## load data
count = read.table("~/03.project/03.circRNA/CircRNAVersion3/02.data/03.workingdir/BtJ.count.xls",
                   header = T,
                   sep = "\t")
head(count)
count = data.frame(row.names = count[,1],
                   count[,-1])
## colData of DEseq
Sample = colnames(count)
Cultivar = gsub(pattern = ".$",replacement = "",Sample)
condition <- factor(Cultivar,levels = unique(Cultivar))
colData = data.frame(row.names = Sample,
                     condition = condition,
                     lib_size = colSums(count),
                     rep = rep(c(1,2,3),times = 9))
## function by Dr. Hu get mid-parent-value
getMid=function(count, coldata, index.a,index.b){
  meanLibSize = mean(coldata$lib_size[c(index.a,index.b)])
  g = expand.grid(index.a,index.b)
  i=1
  mid=(count[,g$Var1[i]]/coldata$lib_size[g$Var1[i]] + count[,g$Var2[i]]/coldata$lib_size[g$Var2[i]])*meanLibSize/2
  for(i in 2:nrow(g)){
    mid=cbind(mid,(count[,g$Var1[i]]/coldata$lib_size[g$Var1[i]] + count[,g$Var2[i]]/coldata$lib_size[g$Var2[i]])*meanLibSize/2 )
  }
  mid=data.frame(mid)
  names(mid)=paste0("mid.",paste0("S",coldata$rep[g$Var1],"-",coldata$rep[g$Var2]))
  return(mid)
}

## parent index
# fiber
BtF = which(colData$condition == "Bt_F")
JF = which(colData$condition == "J_F")
# leaf
BtL = which(colData$condition == "Bt_L")
JL = which(colData$condition == "J_L")
# ovule
BtO = which(colData$condition == "Bt_O")
JO = which(colData$condition == "J_O")
## mid parent
p1 = data.frame(n1 = BtF,n2 = BtL,n3 = BtO)
p2 = data.frame(n1 = JF,n2 = JL,n3 = JO)
mid = list()
for (i in 1:3) {
  mid[[i]] = getMid(count = count, coldata = colData, index.a = p1[,i],index.b = p2[,i])
}

## F1 index
## fiber
BtJF = which(colData$condition == "BtJ_F")
## leaf
BtJL = which(colData$condition == "BtJ_L")
## ovule
BtJO = which(colData$condition == "BtJ_O")
##===================Count Matrix list =================##
P1 = list(count[,BtF],count[,BtL],count[,BtO])
P2 = list(count[,JF],count[,JL],count[,JO])
F1 = list(count[,BtJF],count[,BtJL],count[,BtJO]) 
name = c("BtxJ-F","BtxJ-L","BtxJ-O")
head(F1)
#####================== F1 vs MPV ====================####
## source 
system("mkdir 02.MPV")
for (i in 1:length(name)) { 
  count = cbind(F1[[i]],mid[[i]])
  info = data.frame(sample=names(count), cultivar = rep(c("F1","mid"),times=c(3,9)))
  dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ cultivar)
  res = results(DESeq(dds),contrast=c("cultivar","F1","mid"))
  print( summary(res,alpha=.05) ) 
  res = data.frame(circID = rownames(res),
                   res)
  write.table(res, file=paste("02.MPV/",name[i],"_F1vsMPV.xls",sep = ""),  
              sep="\t",
              quote = F,
              row.names = F)
}

#####================== P1 vs P2 ====================####
system("mkdir 03.P1vsP2")

for (i in 1:length(name)) {
  count = cbind(P1[[i]],P2[[i]])
  info = data.frame(sample=names(count), cultivar = rep(c("P1","P2"),times=c(3,3)))
  dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ cultivar)
  res = results(DESeq(dds),contrast=c("cultivar","P1","P2"))
  print( summary(res,alpha=.05) )
  res = data.frame(circID = rownames(res),
                   res)
  write.table(res, file=paste("03.P1vsP2/",name[i],"_P1vsP2.xls",sep = ""),
              sep="\t",
              quote = F)
}

#####================== F1 vs P1 ====================####
system("mkdir 04.F1vsP1")

for (i in 1:length(name)) {
  count = cbind(F1[[i]],P1[[i]])
  info = data.frame(sample=names(count), cultivar = rep(c("F1","P1"),times=c(3,3)))
  dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ cultivar)
  res = results(DESeq(dds),contrast=c("cultivar","F1","P1"))
  print( summary(res,alpha=.05) )
  res = data.frame(circID = rownames(res),
                   res)
  write.table(res, file=paste("04.F1vsP1/",name[i],"_F1vsP1.xls",sep = ""),
              sep="\t",
              quote = F)
}

#####================== F1 vs P2 ====================####
system("mkdir 05.F1vsP2")

for (i in 1:length(name)) {
  count = cbind(F1[[i]],P2[[i]])
  info = data.frame(sample=names(count), cultivar = rep(c("F1","P2"),times=c(3,3)))
  dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ cultivar)
  res = results(DESeq(dds),contrast=c("cultivar","F1","P2"))
  print( summary(res,alpha=.05) )
  res = data.frame(circID = rownames(res),
                   res)
  write.table(res, file=paste("05.F1vsP2/",name[i],"_F1vsP2.xls",sep = ""),
              sep="\t",
              quote = F)
}
# ## new method
# dds <- DESeqDataSetFromMatrix( countData = count, colData = colData, design = ~ condition)
# 
# pairwiseDE<-function(dds, contrast,savePath)
# {
#   # DE analysis
#   print(contrast)
#   ddsPW <-dds[,dds$condition %in% contrast]
#   ddsPW$condition<-droplevels(ddsPW$condition)
#   res <- results(DESeq(ddsPW), contrast = c("condition",contrast))  # make sure in order
#   print( summary(res,alpha=.05) ) # print results
#   res = data.frame(circID = rownames(res),
#                    res)
#   write.table(res, file=paste(savePath,paste(contrast, collapse="vs"),".txt", sep=""), 
#               sep="\t",
#               quote = F,
#               row.names = F)
# }
# batch <- rbind(
#   c("Bt_L","J_L"),
#   c("Bt_L","SY_L"),
#   c("Bt_F","J_F"),
#   c("Bt_F","SY_F"),
#   c("Bt_O","J_O"),
#   c("Bt_O","SY_O"),
#   c("BtJ_L","J_L"),
#   c("BtSY_L","SY_L"),
#   c("BtJ_F","J_F"),
#   c("BtSY_F","SY_F"),
#   c("BtJ_O","J_O"),
#   c("BtSY_O","SY_O"),
#   c("BtJ_L","Bt_L"),
#   c("BtSY_L","Bt_L"),
#   c("BtJ_F","Bt_F"),
#   c("BtSY_F","Bt_F"),
#   c("BtJ_O","Bt_O"),
#   c("BtSY_O","Bt_O")
# )
# system("mkdir tmp")
# apply(batch,1,function(x) pairwiseDE(dds,x,savePath = "tmp/"))

###################################
## Expression Dominance Analysis ##
###################################

# Parents - P1, P2 or A, D
# Mid parental - Mid
# Hybrid/allopolyploid - T

## if we set padj as the cutoff most groups has DEGs number lower than 10 so, we set the qvalue < 0.05 as the cutoff

classDominance<-function(TvsMid, TvsP1, TvsP2, P1vsP2, log2fc.threshold=0, reverseTvsP =FALSE, Pnames=c("Maternal","Paternal"))
{
  # Hybrid/polyploid vs Mid parental val
  TvsMid <- data.frame(TvsMid[,c("log2FoldChange", "lfcSE", "pvalue")])
  names(TvsMid) <- c("TvsMid", "TvsMid.SE", "TvsMid.pvalue")
  # Hybrid/polyploid vs Parent 1
  TvsP1 <- data.frame(TvsP1[,c("log2FoldChange", "lfcSE", "pvalue")])
  names(TvsP1) <- c("TvsP1", "TvsP1.SE", "TvsP1.pvalue")
  # Hybrid/polyploid vs Parent 2
  TvsP2 <- data.frame(TvsP2[,c("log2FoldChange", "lfcSE", "pvalue")])
  names(TvsP2) <- c("TvsP2", "TvsP2.SE", "TvsP2.pvalue")
  # Parent 1 vs Parent 2
  P1vsP2 <- data.frame(P1vsP2[,c("log2FoldChange", "lfcSE", "pvalue")])
  names(P1vsP2) <- c("P1vsP2", "P1vsP2.SE", "P1vsP2.pvalue")
  
  if(reverseTvsP==TRUE){
    TvsP1$TvsP1 = -TvsP1$TvsP1
    TvsP2$TvsP2 = -TvsP2$TvsP2
  }
  
  tbl = cbind(TvsMid, TvsP1, TvsP2, P1vsP2)
  
  tbl$TvsMid[is.na(tbl$TvsMid)] =0
  tbl$TvsP1[is.na(tbl$TvsP1)] =0
  tbl$TvsP2[is.na(tbl$TvsP2)] =0
  tbl$P1vsP2[is.na(tbl$P1vsP2)] =0
  
  # judge
  tbl$additivity <- ifelse(abs(tbl$TvsMid)>log2fc.threshold & tbl$TvsMid.pvalue<0.2 & !is.na(tbl$TvsMid.pvalue), "T!=Mid", "T=Mid")
  tbl$TvsP1.reg <- ifelse(abs(tbl$TvsP1)>log2fc.threshold & tbl$TvsP1.pvalue<0.2 & !is.na(tbl$TvsP1.pvalue), "T!=P1", "T=P1")
  tbl$TvsP2.reg <- ifelse(abs(tbl$TvsP2)>log2fc.threshold & tbl$TvsP2.pvalue<0.2 & !is.na(tbl$TvsP2.pvalue), "T!=P2", "T=P2")
  tbl$P1vsP2.reg <- ifelse(abs(tbl$P1vsP2)>log2fc.threshold & tbl$P1vsP2.pvalue<0.2 & !is.na(tbl$P1vsP2.pvalue), ifelse(tbl$P1vsP2>log2fc.threshold & tbl$P1vsP2.pvalue<0.2 & !is.na(tbl$P1vsP2.pvalue), "P1>P2","P1<P2"), "P1=P2")
  
  # together
  tbl$class <- paste(tbl$P1vsP2.reg, tbl$additivity, tbl$TvsP1.reg, tbl$TvsP2.reg, sep=";")
  
  # assign category
  tbl$category = "7.Other"
  tbl$category[grep("T=Mid",tbl$class)] = "1.Additivity" # hybrid=Mid while P1!=P2
  tbl$category[grep("P1=P2;T=Mid",tbl$class)] = "6.Conserved"  # P1=P2=Mid=hybrid
  tbl$category[grep("P1>P2;T!=Mid;T=P1;T!=P2",tbl$class)] = paste0("2.",Pnames[1],"-dominant, higher")
  tbl$category[grep("P1<P2;T!=Mid;T=P1;T!=P2",tbl$class)] = paste0("2.",Pnames[1],"-dominant, lower")
  tbl$category[grep("P1>P2;T!=Mid;T!=P1;T=P2",tbl$class)] = paste0("3.",Pnames[2],"-dominant, lower")
  tbl$category[grep("P1<P2;T!=Mid;T!=P1;T=P2",tbl$class)] = paste0("3.",Pnames[2],"-dominant, higher")
  tbl$category[grepl("T!=Mid;T!=P1;T!=P2",tbl$class) & tbl$TvsP1>0 & tbl$TvsP2>0 & tbl$P1vsP2>0] = paste0("4.Transgressive Up: ",Pnames[1]," higher")
  tbl$category[grepl("T!=Mid;T!=P1;T!=P2",tbl$class) & tbl$TvsP1>0 & tbl$TvsP2>0 & tbl$P1vsP2<0] = paste0("4.Transgressive Up: ",Pnames[2]," higher")
  tbl$category[grepl("T!=Mid;T!=P1;T!=P2",tbl$class) & tbl$TvsP1<0 & tbl$TvsP2<0 & tbl$P1vsP2>0] = paste0("5.Transgressive Down: ",Pnames[1]," higher")
  tbl$category[grepl("T!=Mid;T!=P1;T!=P2",tbl$class) & tbl$TvsP1<0 & tbl$TvsP2<0 & tbl$P1vsP2<0] = paste0("5.Transgressive Down: ",Pnames[2]," higher")
  return(tbl)
}

## import data
# TvsMid, TvsP1, TvsP2, P1vsP2
# T vs Mid
path.raw <- "02.MPV/"
fileNames.raw <- dir(path.raw, pattern = "*.xls") 
filePath.raw <- sapply(fileNames.raw, function(x){ 
  paste(path.raw,x,sep='/')})   
TvsMid <- lapply(filePath.raw, function(x){
  read.table(x, sep = "\t",header = TRUE,quote = "",stringsAsFactors = FALSE)}) 
# T vs P1
path.raw <- "04.F1vsP1/"

fileNames.raw <- dir(path.raw, pattern = "*.xls") 
filePath.raw <- sapply(fileNames.raw, function(x){ 
  paste(path.raw,x,sep='/')})   
TvsP1 <- lapply(filePath.raw, function(x){
  read.table(x, sep = "\t",header = TRUE,quote = "",stringsAsFactors = FALSE)}) 
# T vs P2
path.raw <- "05.F1vsP2/"

fileNames.raw <- dir(path.raw, pattern = "*.xls") 
filePath.raw <- sapply(fileNames.raw, function(x){ 
  paste(path.raw,x,sep='/')})   
TvsP2 <- lapply(filePath.raw, function(x){
  read.table(x, sep = "\t",header = TRUE,quote = "",stringsAsFactors = FALSE)}) 
# P1 vs P2
path.raw <- "03.P1vsP2/"

fileNames.raw <- dir(path.raw, pattern = "*.xls") 
filePath.raw <- sapply(fileNames.raw, function(x){ 
  paste(path.raw,x,sep='/')})   
P1vsP2 <- lapply(filePath.raw, function(x){
  read.table(x, sep = "\t",header = TRUE,quote = "",stringsAsFactors = FALSE)}) 
## name
name = gsub(pattern = "_.*",replacement = "",fileNames.raw)

ExpPattern = list()
system("mkdir 06.ExpPattern")
for (i in 1:length(name)) {
  ExpPattern[[i]] = classDominance(TvsMid = TvsMid[[i]], TvsP1 = TvsP1[[i]], 
                                   TvsP2 = TvsP2[[i]], P1vsP2 = P1vsP2[[i]] ,
                                   log2fc.threshold = 0,reverseTvsP = FALSE,Pnames = c("Maternal","Paternal") )
  ExpPattern[[i]] = data.frame(circID = rownames(ExpPattern[[i]]),
                               ExpPattern[[i]])
  write.table(x =ExpPattern[[i]], 
              file = paste("06.ExpPattern/",name[i],".ExpPattern.xls",sep = ""),
              sep = "\t",
              row.names = F,
              quote = F)
}



head(ExpPattern[[1]])



for (i in 1:length(ExpPattern)) {
  categories[[i]] = data.frame(ExpPattern[[i]]$category)
}

x = c("1.Additivity","2.Maternal-dominant, higher","2.Maternal-dominant, lower",
      "3.Paternal-dominant, higher","3.Paternal-dominant, lower",
      "4.Transgressive Up: Maternal higher","4.Transgressive Up: Paternal higher",
      "5.Transgressive Down: Maternal higher","5.Transgressive Down: Paternal higher",
      "6.Conserved","7.Other")
catecolor = data.frame(category = x,
                       color = c("green","blue","#00CED1",
                                 "#FF1493","red",
                                 "orange","#808000",
                                 "gold","purple",
                                 "#C0C0C0","black")
)

## 4-quadrant chart
FourQuadrantChart = function(data,filename,filepath,catecolor){
  ##======1st step=======##
  ## add color information for each category and make a df with category2color
  tmp1 = data.frame(category = data[order(data$category),ncol(data)] %>% unique(.))
  tmp2 = left_join(tmp1,catecolor, by = "category")
  tmp3 = left_join(data,catecolor,by = "category")
  ## main plot
  p = ggplot(tmp3,aes(x = TvsP1, y = TvsP2, color = category)) + 
    geom_hline(aes(yintercept=0), colour="black", linetype="dashed")+
    geom_vline(aes(xintercept = 0), color = "black", linetype = "dashed")+
    geom_point(size = 0.8)+
    scale_color_manual(values = tmp2$color)+
    ylim(-10,10)+ 
    xlim(-10,10)+
    theme_light()+ 
    xlab(expression(paste(log[2],"(",Hybrid/Maternal,")")))+
    ylab(expression(paste(log[2],"(",Hybrid/Paternal,")")))+
    theme(legend.position = "bottom",axis.title = element_text(size = 16),axis.text = element_text(size = 14))
  ##======2nd step============##
  ## add Tags which do not contain additivity conserved and other type.
  ## make a tag df
  tag <- filter(data, category != "1.Additivity" & category != "6.Conserved" & category != "7.Other")
  ## add label layer
  p1 = p + geom_label_repel(data = tag,
                            mapping = aes(x = TvsP1, y = TvsP2, label = circID),
                            arrow = arrow(length = unit(0.01, "npc"), type = "closed", ends = "first"),
                            size = 3,vjust = -3)
  # theme(legend.position="none")
  # ggsave(filename = paste(filepath,filename,"4-quadrant.tiff",sep = ""),plot = p1,width = 12,height = 12,units = 'cm',dpi = 600) 
  # ggsave(filename = paste(filepath,filename,"4-quadrant.svg",sep = ""),plot = p1,width = 12,height = 12,units = 'cm',dpi = 600)
  # ggsave(filename = paste(filepath,filename,"4-quadrant.pdf",sep = ""),plot = p1,width = 12,height = 12,units = 'cm',dpi = 600)
  return(p1)
}
# ## plot
# for (i in 1:length(name)) {
#   FourQuadrantChart(data = ExpPattern[[i]],filename =name[i] ,filepath = "06.ExpPattern/",catecolor = catecolor)
# }

pF = FourQuadrantChart(data = ExpPattern[[1]],filename =name[1] ,filepath = "06.ExpPattern/",catecolor = catecolor)
pL = FourQuadrantChart(data = ExpPattern[[2]],filename =name[2] ,filepath = "06.ExpPattern/",catecolor = catecolor)
pO = FourQuadrantChart(data = ExpPattern[[3]],filename =name[3] ,filepath = "06.ExpPattern/",catecolor = catecolor)

p_all = ggarrange(pF,pL,pO,
          ncol = 3,labels = c("A","B","C"),
          common.legend = TRUE,
          legend = "bottom")
ggsave(plot = p_all,filename = "06.ExpPattern/finalplot.pdf",width = 18,height =6.5 )


