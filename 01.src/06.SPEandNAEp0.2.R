###############################################
#     Prj: Heterosis-circ
#     Assignment: SPE and NAE analysis
#     Author: Shawn Wang
#     Date: Sep 23, 2020
###############################################
## SPE(single parent expressed) and NAE(none additivity expressed) was thought of one of the factor driving heterosis
## Reference ==> Luo, Z., Qian, J., Chen, S., Li, L., 2019. Dynamic patterns of circular and linear RNAs in maize hybrid and parental lines. Theoretical and Applied Genetics.. doi:10007/s00122-019-03489-9
## CPM <0 in one parent and >1 in other parent
## CPM <0 in hybrid and > 1 in at least in one of its parents
##===============01.prepare========================
setwd("~/03.project/03.circRNA/CircRNAVersion3/03.process/06.SPENAE/")
library(dplyr)
library(ggplot2)
library(stringr)
library(edgeR)
library(pheatmap)
library(reshape2)
options(stringsAsFactors = F)
## import count data and transformed as CPM
rawcount = read.table("~/03.project/03.circRNA/CircRNAVersion3/02.data/03.workingdir/BtJ.count.xls",
                      header = T,
                      sep = "\t")
head(rawcount)
clean_count = select(rawcount,contains("Bt"),contains("J"),-contains("SY"),-contains("Z12"))
head(clean_count)
## filter SPE genes
## 1. split into 3 tissues
fiber = select(clean_count,contains("_F"))
leaf = select(clean_count,contains("_L"))
ovule = select(clean_count,contains("_O"))
cpmlist = list(fiber = data.frame(cpm(fiber)),
               leaf = data.frame(cpm(leaf)),
               ovule = data.frame(cpm(ovule)))
head(cpmlist$fiber)
head(cpmlist[[1]])
## 2. function to judge expressed or not with cutoff cpm < 0 & cpm > 1 in at least two of three samples of each cultivar
AddTags = function(cpm){
  for (i in 1:nrow(cpm)) {
    ## Hybrid
    if (sum(cpm[i,1:3] > 1) > 1) {
      cpm[i,10] <- "Ture"
    }else if (sum(cpm[i,1:3] < ) > 1) {
      cpm[i,10] <- "False"
    }else {
      cpm[i,10] <- "None"
    }
    ## Maternal
    if (sum(cpm[i,4:6] > 1) > 1) {
      cpm[i,11] <- "Ture"
    }else if (sum(cpm[i,4:6] < ) > 1) {
      cpm[i,11] <- "False"
    }else {
      cpm[i,11] <- "None"
    }
    ## Paternal
    if (sum(cpm[i,7:9] > 1) > 1) {
      cpm[i,12] <- "Ture"
    }else if (sum(cpm[i,7:9] < ) > 1) {
      cpm[i,12] <- "False"
    }else {
      cpm[i,12] <- "None"
    }
  }
  return(cpm)
}
## 3. based on the sample tag, divide the circRNAs into SPE and NAE categories.
SPEandNAEfilter = function(cpm){
  ## 01. add ture or false tag
  cpm = AddTags(cpm) 
  colnames(cpm)[(ncol(cpm)-2):ncol(cpm)] = c("H.tag","M.tag","P.tag") # change column names
  ## 02. add tags for each circRNA if it belongs to SPE.
  cpm[,13] = case_when(
    cpm$M.tag == "Ture" & cpm$P.tag == "False" ~ "SPE-Maternal",
    cpm$M.tag == "False" & cpm$P.tag == "Ture" ~ "SPE-Paternal"
  )
  ## 03. same as SPE for NAE
  cpm[,14] = case_when(
    cpm$M.tag == "False" & cpm$P.tag == "False" & cpm$H.tag == "Ture" ~ "NAE-Hybrid",
    cpm$M.tag == "Ture" & cpm$P.tag == "Ture" & cpm$H.tag == "False" ~ "NAE-Parent",
    cpm$M.tag == "False" & cpm$P.tag == "Ture" & cpm$H.tag == "False" ~ "NAE-Parent",
    cpm$M.tag == "Ture" & cpm$P.tag == "False" & cpm$H.tag == "False" ~ "NAE-Parent",
  )
  ## finish work, add rownames and circID tag
  colnames(cpm)[(ncol(cpm)-1):ncol(cpm)] = c("SPE","NAE")
  cpm = data.frame(row.names = rawcount$circID,
                   circID = rawcount$circID,
                   cpm)
  return(cpm)
}
fiber.result = SPEandNAEfilter(cpm = cpmlist$fiber)
leaf.result = SPEandNAEfilter(cpm = cpmlist$leaf)
ovule.result = SPEandNAEfilter(cpm = cpmlist$ovule)
spetype = data.frame(circID = rownames(fiber.result),
                     SPE.F = fiber.result$SPE,
                     ENAE.F = fiber.result$NAE,
                     SPE.L = leaf.result$SPE,
                     ENAE.L = leaf.result$NAE,
                     SPE.O = ovule.result$SPE,
                     ENAE.O = ovule.result$NAE)
head(spetype)
head(fiber.result)
fiber.SPE = filter(fiber.result,SPE != "NA")[,c(1,14)]
fiber.NAE = filter(fiber.result,NAE != "NA")[,c(1,15)]
leaf.SPE = filter(leaf.result,SPE != "NA")[,c(1,14)]
leaf.NAE = filter(leaf.result,NAE != "NA")[,c(1,15)]
ovule.SPE = filter(ovule.result,SPE != "NA")[,c(1,14)]
ovule.NAE = filter(ovule.result,NAE != "NA")[,c(1,15)]
fiber.m = full_join(fiber.SPE,fiber.NAE,"circID")
leaf.m = full_join(leaf.SPE,leaf.NAE,"circID")
ovule.m = full_join(ovule.SPE,ovule.NAE,"circID")
x = list(fiber.m,leaf.m,ovule.m)
y = list()
tpm = read.table("~/03.project/03.circRNA/CircRNAVersion3/02.data/03.workingdir/BtJ.tpmNoRep.xls",
                 header = T,
                 sep = "\t")
tpml = list(tpm[,c(1,2,3,4)],
            tpm[,c(1,5,6,7)],
            tpm[,c(1,8,9,10)])
head(tpm)
i = 1
for (i in 1:3) {
  x[[i]]$NAE = gsub(pattern = "NAE",replacement = "ENAE",x = x[[i]]$NAE)
  y[[i]] = data.frame(circID = x[[i]]$circID,
                      type = paste(x[[i]]$SPE,x[[i]]$NAE,sep = ","))
  a = gsub(pattern = ",NA",replacement = "",x = y[[i]]$type)%>%gsub(pattern = "NA,",replacement = "",.)
  y[[i]]$type = factor(a,levels = c("SPE-Maternal","SPE-Paternal,ENAE-Parent",
                             "SPE-Paternal","SPE-Maternal,ENAE-Parent",
                             "ENAE-Hybrid","ENAE-Parent"))
  y[[i]] = inner_join(y[[i]],tpml[[i]],by = "circID")
  y[[i]] = y[[i]][order(y[[i]]$type),]

}

heatNAE = function(data,name){
  x = data.frame(row.names = data$circID,
                 data[,-c(1,2)])
  x = x[,c(2,1,3)]
  annrow = data.frame(row.names = rownames(x),
                      class = data$type)
  pheatmap(mat = log10(x+1),annotation_row  = annrow,scale = "row",
           show_rownames = F,cutree_cols = 3,
           colorRampPalette(c("green", "black", "red"))(50),border_color = "NA",
           cluster_rows = F,cluster_cols = F,filename = paste("~/03.project/03.circRNA/CircRNAVersion3/03.process/06.SPENAE/",name,".pdf",sep = ""),width = 8, height = 9)
}


# pheatmap(mat = log10(x+1),annotation_row  = annrow,scale = "row",
#          show_rownames = T,cutree_cols = 3,fontsize_row = 5,
#          colorRampPalette(c("green", "black", "red"))(50),border_color = "NA",
#          cluster_rows = F,cluster_cols = T)

heatNAE(y[[1]],"fiber")
heatNAE(y[[2]], "leaf")
heatNAE(y[[3]], "ovule")
write.table(spetype,file = "spetype.xls",
            row.names = F,
            sep = "\t",
            quote = F)
## TPM box-jitter plot
library(tidyverse)
library(cowplot)
library(ggpol)
library(grid)
library(rlang)
cultivar = c("Bt","BtJ","J")
tissue = c("F","L","O")
i  = 3
for (i in 1:length(y)) {
  a = melt(y[[i]][,-1],id.vars = "type")
  levels = paste(cultivar,tissue[i],sep = "_")
  a$variable = factor(x = a$variable, levels = levels)
  p = ggplot(a)+
    geom_boxjitter(
      aes( y = value, fill = type),
      jitter.shape = 21, jitter.color = NA, jitter.size = 1, 
      jitter.params = list(height = 0, width = 0.04),
      outlier.color = NA, errorbar.draw = TRUE
    )+
    facet_wrap(facets = ~variable)+
    scale_fill_manual(values = c("#ff83ff","#ff9289","#82b7ff","#00d65c","#00dae0","#d3ba00"))+
    ylim(... = c(0,7000))+
    theme_bw()+
    theme(title = element_blank(),
          legend.position = "none",
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 14))
  ggsave(p,filename = paste("~/03.project/03.circRNA/CircRNAVersion3/03.process/06.SPENAE/boxjitter_",
                            tissue[i],".pdf",sep = ""),width = 8, height = 5,dpi = 300 )
}













## 4. Merge them in a list.
SPE = list(SPE.M.fiber = data.frame(circID = filter(fiber.result,SPE == "SPE-Maternal")[,1]),
           SPE.P.fiber = data.frame(circID = filter(fiber.result,SPE == "SPE-Paternal")[,1]),
           SPE.M.leaf = data.frame(circID = filter(leaf.result,SPE == "SPE-Maternal")[,1]),
           SPE.P.leaf = data.frame(circID = filter(leaf.result,SPE == "SPE-Paternal")[,1]),
           SPE.M.ovule = data.frame(circID = filter(ovule.result,SPE == "SPE-Maternal")[,1]),
           SPE.P.ovule = data.frame(circID = filter(ovule.result,SPE == "SPE-Paternal")[,1]))
NAE = list(NAE.H.fiber = data.frame(circID = filter(fiber.result,NAE == "NAE-Hybrid")[,1]),
           NAE.P.fiber = data.frame(circID = filter(fiber.result,NAE == "NAE-Parent")[,1]),
           NAE.H.leaf = data.frame(circID = filter(leaf.result,NAE == "NAE-Hybrid")[,1]),
           NAE.P.leaf = data.frame(circID = filter(leaf.result,NAE == "NAE-Parent")[,1]),
           NAE.H.ovule = data.frame(circID = filter(ovule.result,NAE == "NAE-Hybrid")[,1]),
           NAE.P.ovule = data.frame(circID = filter(ovule.result,NAE == "NAE-Parent")[,1]))
circ2host = read.table("~/03.Project/03.circRNA/CircRNAVersion3/02.data/03.workingdir/circ2host",
                       header = T,
                       sep = "\t")
tissue = c("F","F","L","L","O","O")
category1 = rep(c("Maternal","Paternal"),times = 3)
category2 = rep(c("Hybrid","Parent"),times = 3)
SPE.mRNA = list()
for (i in 1:length(SPE)) {
  SPE.mRNA[[i]] = inner_join(SPE[[i]],circ2host,by = "circID")
  out = data.frame(geneID = SPE.mRNA[[i]][,2])
  write.table(x = out,file = paste("~/03.Project/03.circRNA/CircRNAVersion3/03.process/06.SPENAE/SPE.",tissue[i],"_",category1[i],".xls",sep = ""),
              row.names = F,
              quote = F,
              sep = "\t")
}
NAE.mRNA = list()
for (i in 1:length(NAE)) {
  NAE.mRNA[[i]] = inner_join(NAE[[i]],circ2host,by = "circID")
  out = data.frame(geneID = NAE.mRNA[[i]][,2])
  write.table(x = out,file = paste("~/03.Project/03.circRNA/CircRNAVersion3/03.process/06.SPENAE/NAE.",tissue[i],"_",category2[i],".xls",sep = ""),
              row.names = F,
              quote = F,
              sep = "\t")
}
save(NAE.mRNA,SPE.mRNA,file = "~/03.Project/03.circRNA/CircRNAVersion3/03.process/06.SPENAE/SPENAEmap.Rdata")

## It is hard to get enrichment result because the host gene set is too small.
## 5. out put and do go analysis with TBtools
## GO enrichment result were trimmed by REViGO
## import trimmed GO result
# path.raw <- "~/03.Project/03.circRNA/CircRNAVersion3/03.process/06.SPENAE/result/"
# fileNames.raw <- dir(path.raw, pattern = "*.txt") 
# filePath.raw <- sapply(fileNames.raw, function(x){ 
#   paste(path.raw,x,sep='/')})   
# goTrimmed <- lapply(filePath.raw, function(x){
#   read.table(x, sep = "\t",header = TRUE,quote = "",stringsAsFactors = FALSE)}) 
# head(goTrimmed[[1]])
# ## 6. dotplot ==> top ranked 15 terms of each group
# name <- gsub(pattern = ".GO.txt",replacement = "",x = fileNames.raw)
# name
# colnames(goTrimmed[[1]])
# godotplot = function(left,right,lname,rname,name){
#   ## top 15 terms
#   x = left[order(left$P_value)[1:15],]
#   y = right[order(right$P_value)[1:15],]
#   ## Gene ratio calculate
#   x$GeneRatio = x$HitsGenesCountsInSelectedSet/x$AllGenesCountsInSelectedSet
#   y$GeneRatio = y$HitsGenesCountsInSelectedSet/y$AllGenesCountsInSelectedSet
#   ## ordered the GO name for plot
#   x = x[order(x$GeneRatio),]
#   y = y[order(y$GeneRatio,decreasing = T),]
#   ## data cleanning
#   leftdata = data.frame(GO_Name = factor(x$GO_Name,levels = x$GO_Name),
#                         P_value = -log10(x$P_value),
#                         GeneRatio = x$GeneRatio,
#                         Cultivar = rname)
#   rightdata = data.frame(GO_Name = factor(y$GO_Name,levels = y$GO_Name),
#                          P_value = -log10(y$P_value),
#                          GeneRatio = y$GeneRatio,
#                          Cultivar = lname)
#   ## merge
#   z = rbind(leftdata,rightdata)
#   ## plot
#   p = ggplot(data = z, mapping = aes(x = Cultivar,y = GO_Name))+
#     geom_point(aes(size = GeneRatio,color = P_value))+
#     scale_color_gradient(low = "blue", high = "red", guide=guide_colorbar(reverse=TRUE))+
#     scale_y_discrete(labels=function(y) stringr::str_wrap(y,width=40)) +
#     theme_bw() +
#     theme(axis.text.x = element_text(colour = "black",
#                                      size = 14, vjust =1 ),
#           axis.text.y = element_text(colour = "black",
#                                      size = 14, hjust =1 ),
#           axis.title = element_text(margin=margin(10, 5, 0, 0),
#                                     color = "black",size = 16),
#           axis.title.y = element_text(angle=90))+ 
#     labs(color=expression(-log[10](Qvalue)),
#          size="Gene Ratio",
#          x = "",
#          y = "")
#   ggsave(p,filename = paste("~/02.Project/03.circRNA/CircRNAVersion2/03.process/06.SPE/02.enrichment/",name,".pdf",sep = ""),width = 9,height = 10,dpi = 300)
# }
# 
# ## plot
# godotplot(left = goTrimmed[[1]],
#           right = goTrimmed[[2]],
#           lname = "Hybrid",
#           rname = "MPV",
#           name = "NAE_F")
# godotplot(left = goTrimmed[[3]],
#           right = goTrimmed[[4]],
#           lname = "Hybrid",
#           rname = "MPV",
#           name = "NAE_L")
# godotplot(left = goTrimmed[[5]],
#           right = goTrimmed[[6]],
#           lname = "Hybrid",
#           rname = "MPV",
#           name = "NAE_O")
# godotplot(left = goTrimmed[[7]],
#           right = goTrimmed[[8]],
#           lname = "Maternal",
#           rname = "Paternal",
#           name = "SPE_F")
# godotplot(left = goTrimmed[[9]],
#           right = goTrimmed[[10]],
#           lname = "Maternal",
#           rname = "Paternal",
#           name = "SPE_L")
# godotplot(left = goTrimmed[[11]],
#           right = goTrimmed[[12]],
#           lname = "Maternal",
#           rname = "Paternal",
#           name = "SPE_O")
# ## plot for the figure 3 main plot
# 
# SPEALL = list(fiber = rbind(SPE.mRNA[[1]],SPE.mRNA[[2]]),
#               leaf = rbind(SPE.mRNA[[3]],SPE.mRNA[[4]]),
#               ovule = rbind(SPE.mRNA[[5]],SPE.mRNA[[6]]))
# NAEALL =  list(fiber = rbind(NAE.mRNA[[1]],NAE.mRNA[[2]]),
#                leaf = rbind(NAE.mRNA[[3]],NAE.mRNA[[4]]),
#                ovule = rbind(NAE.mRNA[[5]],NAE.mRNA[[6]]))


