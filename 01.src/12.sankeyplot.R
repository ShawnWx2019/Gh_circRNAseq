###############################################
#     Prj: Heterosis-circ
#     Assignment: Sankey plot for ceRNA
#     Author: Shawn Wang
#     Date: Dec 24, 2020
###############################################
## for sankey plot data cleanning
setwd("/Volumes/Samsung_T5/03.circRNA/CircRNAVersion3/03.process/10.SankeyPlot/")
library(dplyr)
library(ggplot2)
library(stringr)
library(networkD3)
options(stringsAsFactors = F)
source("~/02.MyScript/BioScript/01.R/ShawnRToolkit.R")
#####=======01. import data=========
## 1 import gene2goid
path = "/Volumes/Samsung_T5/03.circRNA/CircRNAVersion3/03.process/10.SankeyPlot/"
g2go = BatchReadTable(path = path,pattern = "01.*",sep = "\t",header = F,quote = NULL,stringsAsFactors = F)
## name
name = gsub(pattern = "Gene2GOterm.txt",replacement = "",g2go$Name) %>%
  gsub(pattern = "01.",replacement = "",.)
## table
g2gotable = g2go$file
#head(g2gotable)
nodecolor = c("#FF69B4","#00FA9A","#FFD700")
for (i in 1:3) {
  names(g2gotable[[i]]) = c("GID","GO_Name","Tissue")
  g2gotable[[i]]$color = nodecolor[i]
  
}

## import gene2term
g2name = BatchReadTable(path = path, pattern = "02.", sep = "\t", header = T, quote = NULL, stringsAsFactors = F)
## replace the NA as GID
g2nameClean = function(data){
  data[is.na(data)] <- 'aaaa'
  data = data[,-1]
  names(data) = c("GID","G_Name")
  for (i in 1:nrow(data)) {
    if (data[i,2] == 'aaaa') {
      data[i,2] = data[i,1]
    }
  }
  return(data)
}
g2nameTable = list()
tissue = c("fiber","leaf","ovule")
for (i in 1:length(g2name$file)) {
  g2nameTable[[i]] = g2nameClean(g2name$file[[i]])
}
raw.relations <- read.table("/Volumes/Samsung_T5/03.circRNA/CircRNAVersion3/03.process/09.finalpathway/03.other/expmi2mRNA.xls",
                            header = T,
                            sep = "\t")
#head(raw.relations)
## relation split
relation = list(
  fiber = dplyr::filter(raw.relations,Tissue == "Fiber"),
  leaf = dplyr::filter(raw.relations,Tissue == "Leaf"),
  ovule = dplyr::filter(raw.relations,Tissue == "ovule")
)

for (i in 1:3) {
  relation[[i]]$color = nodecolor[i]
}
merge = list()
for (i in 1:3) {
  merge[[i]] = dplyr::inner_join(g2nameTable[[i]], relation[[i]], by = "GID")%>%
    inner_join(g2gotable[[i]][,c(1,2)],by = 'GID')
}
all = bind_rows(merge)
link.T2circ = data.frame(unique(all[,c("Tissue","circID")]),
                         score = 10)
names(link.T2circ) = c("source","target","score")
link.c2mi = data.frame(unique(all[,c("circID","miID")]),
                       score = 10)
names(link.c2mi) = c("source","target","score")
link.mRNA2Gname = data.frame(unique(all[,c("GID","G_Name")]),
                             score = 10)
names(link.mRNA2Gname) = c("source","target","score")
link.mi2G = unique(all[,c("miID","GID")])
names(link.mi2G) = c("source","target")
link.name2go = data.frame(unique(all[,c("G_Name","GO_Name")]),
                          score = 10)
names(link.name2go) = c("source","target","score")
sankeyNode = function(data,tag){
  link.circ2miRNA = data.frame(source = data$circID,
                               target = data$miID,
                               score = 6)
  link.miRNA2GID = data.frame(source = data$miID,
                              target = data$GID,
                              score = 6)
  link.GID2Gname = data.frame(source = data$GID,
                              target = data$G_Name,
                              score = 10)
  link.GName2GOName = data.frame(source = data$G_Name,
                                 target = data$GO_Name,
                                 score = 10)
  link = rbind(link.circ2miRNA,
               link.miRNA2GID,
               link.GID2Gname,
               link.GName2GOName)
  link = unique(link)
  node1 = unique(link$source)
  node2 = unique(link$target)
  node = c(node1,node2)
  node = unique(node)
  node.all = data.frame(node = node)
  write.table(node.all,file = paste0(tag,".xls"),
              row.names = F,
              quote = F,
              sep = "\t")
  x = list(node = node,
           link = link)
  return(x)
}

sankeyPlot = function(path,node,link,height){
  node.all = read.table(path,
                         sep = "\t",
                         header = T)
  node.source = data.frame(source = node,
                           sourceID = seq(1:length(node))-1)
  node.target = data.frame(target = node,
                           targetID = seq(1:length(node))-1)
  
  link = left_join(link,node.source,by = "source") %>%
    left_join(.,node.target,by = "target")
  head(link)
  
  p = sankeyNetwork(Links = link, Nodes = node.all, Value = "score",
                Source = "sourceID", Target = 'targetID',NodeGroup = "group",
                NodeID = "node",nodeWidth = 20, fontSize = 12,height = height, width = 1000)
  return(p)
}

# fiber = sankeyNode(data = merge[[1]],tag = "fibernode")
# leaf = sankeyNode(data = merge[[2]],tag = "leafnode")
# ovule = sankeyNode(data = merge[[3]],tag = "ovulenode")
# p.fiber = sankeyPlot(path = "fibernode.xls",
#                      node = fiber$node,
#                      link = fiber$link,
#                      height = 800 )
## The picture size is too large to fit in the paper， So I want to remove the information of GeneID and Gene Name. Change the score value as the edge number
## node type:
# 01. ==> tisue: fiber leaf ovule
# 02. ==> key circRNAs: with miRNA binding site, without miRNA binding site
# 03. ==> miRNAs: 
# 04. ==> miRNA-mRNA targets: key GO terms, other GO terms, with out GO terms.
# 05. ==> GO terms.

link.T2circ

link.T2circ$lgroup = rep(c("Fiber","Leaf","Ovule"),times = c(2,1,5))
link.T2circ$source = link.T2circ$lgroup
link.T2circ$score = 1
link.c2mi$lgroup = link.c2mi$source
link.c2mi$score = 1
link.mi2G
# count the number of miRNA targets.
tmp1 = plyr::count(link.mi2G,vars = "source")
link.N.mi2G = data.frame(source = tmp1$source,
                       target = paste0(tmp1$source,"-Targets"),
                       score = tmp1$freq,
                       lgroup = tmp1$source)
link.g2GO = data.frame(all[,c("GID","GO_Name")])
#head(link.g2GO)
names(link.g2GO) = c("target",'GO_Name')
tmp2 = left_join(link.mi2G,link.N.mi2G,by = "source")
tmp3 = unique(tmp2[,c("target.x","target.y")])
names(tmp3) = c("target","targetNew")
tmp4 = inner_join(link.g2GO,tmp3,by = "target")
tmp4 = data.frame(source = tmp4$targetNew,
                  target = tmp4$GO_Name)
link.N.G2GO = plyr::count(tmp4)
link.N.G2GO$lgroup = link.N.G2GO$target
names(link.N.G2GO) = c("source","target","score","lgroup")

final = rbind(link.T2circ,link.c2mi,link.N.mi2G,link.N.G2GO)
View(final)
getsankeydata = function(link.all){
  node1 = unique(link.all$source)
  node2 = unique(link.all$target)
  node = c(node1,node2)
  node = unique(node)
  node.all = data.frame(node = node)
  node.source = data.frame(source = node,
                           sourceID = seq(1:length(node))-1)
  node.target = data.frame(target = node,
                           targetID = seq(1:length(node))-1)
  
  linkall = left_join(link.all,node.source,by = "source") %>%
    left_join(.,node.target,by = "target")
  x = list(node = node.all,
           link = linkall)
  return(x)
}
data = getsankeydata(link.all = final)

p = sankeyNetwork(Links = data$link, Nodes = data$node,Value = "score",LinkGroup = "lgroup",
                  Source = "sourceID", Target = 'targetID',
                  NodeID = "node",nodeWidth = 20, fontSize = 12,height = 600, width = 1200)
saveNetwork(p, "./test.html", selfcontained = TRUE)

