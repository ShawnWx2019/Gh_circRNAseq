####################################################
#       Prj: circRNA
#       Assginment: relationships
#       Author: Shawn Wang
#       Date: Dec 12, 2020
####################################################

###############################################################
### circRNAs which was belongs to each expression pattern   ###
### was marked as heterosis-related highlight,              ###
### the network was constructed by the                      ###
### relationships among circRNA - target miRNAs-target mRNA ###
### Todo list:                                              ###
### 01. miRNA 2 circRNA                                     ###
### 02. miRNA 2 mRNA                                        ###
### 03. ppi based on the GO enrichment analysis             ###
###############################################################
##========== 01. import packages =========================
options(stringsAsFactors = F)
library(dplyr)
library(stringr)
library(ggplot2)
##========== 02. relationships ===========================
setwd("/Volumes/Samsung_T5/03.circRNA/CircRNAVersion3/03.process/09.finalpathway/")
raw.relations <- read.table("/Volumes/Samsung_T5/03.circRNA/CircRNAVersion3/03.process/09.finalpathway/03.other/expmi2mRNA.xls",
                            header = T,
                            sep = "\t")
raw.relations[1:4,1:4]
unique(raw.relations[,c(1,2,3)])
fiber = filter(raw.relations,Tissue == "Fiber")
leaf = filter(raw.relations,Tissue == "Leaf")
ovule = filter(raw.relations,Tissue == "ovule")
##============IDmapping for string==============
## 01.fiber
unique(fiber[,c(1,3)])
fiber.mi2mRNA = unique(fiber[,c(3,4)])
## mi2genename
STmap = read.delim("/Volumes/Samsung_T5/03.circRNA/CircRNAVersion3/03.process/09.finalpathway/02.ceRNA/03.stringoutput/01.fiber/fiberMapping.tsv",
                   header = T,
                   sep  = "\t",
                   quote = "")
head(STmap)
STmap.clean = STmap[,c(2,3,6)]
head(STmap.clean)
names(STmap.clean)[1] = "GID"
fiber.mi2NAME = inner_join(fiber.mi2mRNA,STmap.clean,by = "GID")
head(fiber.mi2NAME)
fiber.mi2NAME.out = unique(fiber.mi2NAME[,c(1,4)])
write.table(x = fiber.mi2NAME.out,file = "/Volumes/Samsung_T5/03.circRNA/CircRNAVersion3/03.process/09.finalpathway/02.ceRNA/03.stringoutput/01.fiber/mi2name.tsv",
            row.names = F,
            quote = F,
            sep = "\t")
## 02.ovule
ovule.circ2mi = unique(ovule[,c(1,3)])
ovule.mi2mRNA = unique(ovule[,c(3,4)])
## mi2genename
write.table(ovule.circ2mi,file = "/Volumes/Samsung_T5/03.circRNA/CircRNAVersion3/03.process/09.finalpathway/02.ceRNA/03.stringoutput/03.ovule/circ2mi.tsv",
            row.names = F,
            quote = F,
            sep = "\t")
STmap = read.delim("/Volumes/Samsung_T5/03.circRNA/CircRNAVersion3/03.process/09.finalpathway/02.ceRNA/03.stringoutput/03.ovule/ovule.map.tsv",
                   header = T,
                   sep  = "\t",
                   quote = "")
head(STmap)
STmap.clean = STmap[,c(2,3,6)]
head(STmap.clean)
names(STmap.clean)[1] = "GID"
ovule.mi2NAME = inner_join(ovule.mi2mRNA,STmap.clean,by = "GID")
head(ovule.mi2NAME)
ovule.mi2NAME.out = unique(ovule.mi2NAME[,c(1,4)])

write.table(x = ovule.mi2NAME.out,file = "/Volumes/Samsung_T5/03.circRNA/CircRNAVersion3/03.process/09.finalpathway/02.ceRNA/03.stringoutput/03.ovule/mi2name.tsv",
            row.names = F,
            quote = F,
            sep = "\t")
## 03.leaf
leaf.circ2mi = unique(leaf[,c(1,3)])
leaf.mi2mRNA = unique(leaf[,c(3,4)])
## mi2genename
write.table(leaf.circ2mi,file = "/Volumes/Samsung_T5/03.circRNA/CircRNAVersion3/03.process/09.finalpathway/02.ceRNA/03.stringoutput/02.leaf/circ2mi.tsv",
            row.names = F,
            quote = F,
            sep = "\t")
STmap = read.delim("/Volumes/Samsung_T5/03.circRNA/CircRNAVersion3/03.process/09.finalpathway/02.ceRNA/03.stringoutput/02.leaf/leafmapping.tsv",
                   header = T,
                   sep  = "\t",
                   quote = "")
head(STmap)
STmap.clean = STmap[,c(2,3,6)]
head(STmap.clean)
names(STmap.clean)[1] = "GID"
leaf.mi2NAME = inner_join(leaf.mi2mRNA,STmap.clean,by = "GID")
head(leaf.mi2NAME)
leaf.mi2NAME.out = unique(leaf.mi2NAME[,c(1,4)])

write.table(x = leaf.mi2NAME.out,file = "/Volumes/Samsung_T5/03.circRNA/CircRNAVersion3/03.process/09.finalpathway/02.ceRNA/03.stringoutput/02.leaf/mi2name.tsv",
            row.names = F,
            quote = F,
            sep = "\t")
##############==end==########################
exportxls = function(x,path,name){
  write.table(x = x,
              file = paste0(path,name,".xls"),
              row.names = F,
              quote = F,
              sep = "\t")
}
path = "/Volumes/Samsung_T5/03.circRNA/CircRNAVersion3/03.process/09.finalpathway/"
exportxls(x = data.frame(GID = fiber$GID), path = path, name = "FibermRNAID")
exportxls(x = data.frame(GID = leaf$GID), path = path, name = "LeafmRNAID")
exportxls(x = data.frame(GID = ovule$GID), path = path, name = "OvulemRNAID")
# system("mkdir 01.enrichment")
unique(ovule[,c(1,3)])
## merge the table
BatchReadTable = function(path, pattern,sep = "\t", header = TRUE){
  fileNames.raw <- dir(path, pattern = pattern) 
  filePath.raw <- sapply(fileNames.raw, function(x){ 
    paste(path,x,sep='/')})   
  data.raw <- lapply(filePath.raw, function(x){
    read.table(x, sep = sep,header = header,quote = "",stringsAsFactors = FALSE)}) 
  x = list(Name = fileNames.raw,
           file = data.raw)
  return(x)
}
raw = BatchReadTable(path = "/Volumes/Samsung_T5/03.circRNA/CircRNAVersion3/03.process/09.finalpathway/01.enrichment/01.GOEnrich/01.final/",
                     pattern = "*final.xls")
raw.Table = raw$file
revigo = BatchReadTable(path = "/Volumes/Samsung_T5/03.circRNA/CircRNAVersion3/03.process/09.finalpathway/01.enrichment/03.REViGO/",
                        pattern = "*all.txt*")

revigo.table = revigo$file
# head(revigo.table$fiber.all.txt)
name = revigo$Name
name = gsub(pattern = ".all.txt",replacement = "",x = name)
# View(raw.Table$Fiber.GO.Enrichment.final.xls)
merge <- list()
system("mkdir 04.FinalEnrichTable")
for (i in 1:3) {
  ## change col id
  a1 = revigo.table[[i]]
  b1 = raw.Table[[i]]
  a1 = data.frame(GO_ID = a1$term_ID,
                 plot_X = a1$plot_X,
                 plot_Y = a1$plot_Y)
  merge[[i]] = left_join(a1,b1,by = "GO_ID")
  write.table(x = merge[[i]],file = paste0("/Volumes/Samsung_T5/03.circRNA/CircRNAVersion3/03.process/09.finalpathway/04.FinalEnrichTable/",name[i],".xls"),
              sep = "\t",
              quote = F,
              row.names = F)
}


