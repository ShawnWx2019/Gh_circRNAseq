############################################
#     Prj: CircRNA
#     Assignment: Different type of circRNA
#     Date: 29 Oct, 2020
#     Author: Shawn Wang
############################################
options(stringsAsFactors = F)
library(stringr)
library(dplyr)
library(RColorBrewer)
## 01. table
setwd("~/03.project/03.circRNA/CircRNAVersion3/03.process/01.Indentification/")
raw = read.delim("~/03.Project/03.circRNA/CircRNAVersion3/02.data/03.workingdir/BtJ_circRNAs.txt",
                 header = T,
                 sep = "\t")
head(raw)
clean = raw[,c(1:4,6,7,8,9)]
head(clean)
exonic = filter(clean,str_detect(string = clean$feature, "exon"))
## intronic and intergenic
other = filter(clean,str_detect(string = clean$feature, "in"))
exonic$ExonNum = str_count(string = exonic$feature,pattern = "exon")
other$type = gsub(pattern = ":.*",replacement = "",x = other$feature)
head(other)
exonsum = case_when(
  exonic$ExonNum == 1 ~ "1 exon",
  exonic$ExonNum == 2 ~ "2 exon",
  exonic$ExonNum == 3 ~ "3 exon",
  exonic$ExonNum > 3 ~ "> 3 exon",
  TRUE ~ as.character(exonic$ExonNum)
)
othersum = data.frame(type = other$type)
exonsum = data.frame(type = exonsum)
sum = rbind(othersum,exonsum)
piedata = unique(add_count(sum,type))
head(piedata)
piedata$type = factor(piedata$type,levels = piedata$type)
## plot
label = paste(piedata$type,"(",round(piedata$n/sum(piedata$n)*100),"%)",
              sep = "") 
pdf(file = "~/03.Project/03.circRNA/CircRNAVersion3/04.result/01.Indentification/01.Pie.pdf",
    width = 12,height = 7)
pie(x = piedata$n,labels = label, col = brewer.pal(6,"Set3"),border = "white", cex = 3,radius = 1)
dev.off()
