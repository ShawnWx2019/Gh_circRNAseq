###########################################
#     Prj: circRNA
#     Assignment: Pvalue selection
#     Date: May 14, 2021
#     Author: Shawn Wang
###########################################
setwd("/Volumes/Samsung_T5/03.circRNA/CircRNAVersion3/03.process/02.ExpressionPattern/07.pvalueTest/")
library(tidyverse)
library("ggprism")
library(patchwork)
## The relationships of DECs number and pvalue


## import DESeq2 result

source("~/02.MyScript/BioScript/01.R/ShawnRToolkit.R")
file = BatchReadTable(path = "./",pattern = ".xls",sep = "\t",header = T,quote = NULL,stringsAsFactors = F)
CompTag = gsub(pattern = ".xls",replacement = "",file$Name)
DECres = file$file
DECres




phist <- function(DEresult,tissue,group){
  myData = DEresult$pvalue
  
  gg_b <- ggplot_build(
    ggplot() + geom_histogram(aes(x = myData), binwidth=.1)
  )
  
  result = gg_b$data[[1]]
  p = ggplot(result,mapping = aes(x = x,y = count)) + 
    geom_col(mapping = aes(fill = count))+
    scale_fill_gradient(low = "blue",high = "salmon")+
    theme_bw()+
    xlab("P value")+
    ylab("Number")+
    theme_prism(border = T,base_size = 9)+
    ggtitle(paste(tissue,group,sep  = " - "))
  return(p)
}

p1 = phist(DEresult = DECres$`BtxJ-F_F1vsP1.xls`,tissue = "Fiber",group = "F1vsP1")
p2 = phist(DEresult = DECres$`BtxJ-F_F1vsP2.xls`,tissue = "Fiber",group = "F1vsP2")
p3 = phist(DEresult = DECres$`BtxJ-F_P1vsP2.xls`,tissue = "Fiber",group = "P1vsP2")
p4 = phist(DEresult = DECres$`BtxJ-L_F1vsP1.xls`,tissue = "Leaf",group = "F1vsP1")
p5 = phist(DEresult = DECres$`BtxJ-L_F1vsP2.xls`,tissue = "Leaf",group = "F1vsP2")
p6 = phist(DEresult = DECres$`BtxJ-L_P1vsP2.xls`,tissue = "Leaf",group = "P1vsP2")
p7 = phist(DEresult = DECres$`BtxJ-O_F1vsP1.xls`,tissue = "Ovule",group = "F1vsP1")
p8 = phist(DEresult = DECres$`BtxJ-O_F1vsP2.xls`,tissue = "Ovule",group = "F1vsP2")
p9 = phist(DEresult = DECres$`BtxJ-O_P1vsP2.xls`,tissue = "Ovule",group = "P1vsP2")

p = (p1+p2+p3)/(p4+p5+p6)/(p7+p8+p9)+
  plot_annotation(
    title = "P value distribution"
  )+plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
ggsave(p,filename = "/Volumes/Samsung_T5/杂优投稿文件/03.MajorRevision/Pvalue-distribution.pdf",width = 10,height = 7)



PvalPltData = function(res,tissue,group,comp){
  pvec = seq(0,1,0.001)
  num = 0
  for (i in 1:length(pvec)) {
    res %>% 
      filter(pvalue <= pvec[i]) %>% 
      nrow(.) -> num[i]
  }
  
  plt.data = data.frame(
    pvalue = pvec,
    DECnum = num,
    tissue = tissue,
    Compare = comp,
    group = group
  )
  return(plt.data)
}

CompTag
Tissue = rep(c("Fiber","Leaf","Ovule"),each = 3)
CompType = rep(c("F1 vs P1","F1 vs P2", "P1 vs P2"),times = 3)
plt.list = list()
for (i in 1:length(CompTag)) {
  plt.list[[i]] = PvalPltData(res = DECres[[i]],
                              tissue = Tissue[i],
                              comp = CompTag[i],
                              group = CompType[i])
}

plt.data = bind_rows(plt.list)
head(plt.data)

densityPlot = function(plt,tissue){
  ggplot(plt,aes(x = pvalue,y = DECnum))+
    geom_line(aes(color = group))+
    scale_color_manual(values = c("salmon","navy","gold"))+
    geom_vline(xintercept = 0.05,color = "blue",size = 0.8,linetype = "dashed",alpha = 2)+
    geom_vline(xintercept = 0.2,color = "red",size = 0.8,linetype = "dashed",alpha = 2)+
    theme_bw()+
    theme_prism(
      border = T
    )+
    geom_hline(yintercept = max(plt$DECnum)*0.1,color = "green",size = 0.8,linetype = "dashed",alpha = 2)+
    scale_x_continuous(
      limits = c(0,1),
      breaks = seq(0,1,0.2)
    )+
    theme(
      legend.position = "bottom"
    )+ 
    ylim(c(0,340))+
    ggtitle(tissue)+
    theme(axis.title = element_text(size = 12),
          plot.title = element_text(size = 13))+
    xlab(label = "P-value")+
    ylab(label = "DEC number")
}

fiber = plt.data %>% 
  filter(tissue == "Fiber") %>% 
  densityPlot(plt = .,tissue = "Fiber")


leaf = plt.data %>% 
  filter(tissue == "Leaf") %>% 
  densityPlot(plt = .,tissue = "Leaf")

ovule = plt.data %>% 
  filter(tissue == "Ovule") %>% 
  densityPlot(plt = .,tissue = "Ovule")


p = (fiber+leaf+ovule) +
  plot_annotation(
    title = "Relationships between DEC numbers and P-value",
    caption = "The blue vertical line represents P-value = 0.05.\n The red vertical line represents P-value = 0.2.\n The green horizontal line represents 10% of the total number of expressed circRNA in the corresponding tissue"
      )+plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
ggsave(p,filename = "Pvalue-DECnum.pdf",width = 10,height = 7)
