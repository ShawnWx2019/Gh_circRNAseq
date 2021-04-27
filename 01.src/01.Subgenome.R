###########################################
#       Prj: Heterosis ==> circRNA
#       Assignment: Tissue specific & subgenome all
#       Author: Shawn Wang
#       Date: Sep 11,2020
###########################################
## get circRNA with non-zero read count of each tissue and labeled as tissue detected circRNAs
library(ggplot2)
library(dplyr)
library(stringr)
options(stringsAsFactors = F)
tpm = read.table("~/03.project/03.circRNA/CircRNAVersion3/02.data/03.workingdir/BtJ.tpm.xls",
                 header = T,
                 sep = "\t")
head(tpm)
fiber = data.frame(circID = tpm$circID,
                   select(tpm,contains("_F")))
leaf = data.frame(circID = tpm$circID,
                  select(tpm,contains("_L")))
ovule = data.frame(circID = tpm$circID,
                   select(tpm,contains("_O")))
## remove zero count circRNAs of each tissue 
fiber1 = fiber[apply(fiber[,-1], 1, sum) != 0,]
leaf1 = fiber[apply(leaf[,-1], 1, sum) != 0,]
ovule1 = ovule[apply(ovule[,-1], 1, sum) != 0,]
## subgenome data
info = read.delim("~/03.project/03.circRNA/CircRNAVersion3/02.data/02.CIRIresult/seq_infor/novel_circRNAs.txt",
                  header = T,
                  sep = "\t")
info.clean = info[,c(1,2)]
names(info.clean) = c("circID","chr")
tmp.F = left_join(data.frame(circID = fiber1[,1]),info.clean,by = "circID")
tmp.L = left_join(data.frame(circID = leaf1[,1]),info.clean,by = "circID")
tmp.O = left_join(data.frame(circID = ovule1[,1]),info.clean,by = "circID")
## count sub chr number
#Tiss = c("fiber","leaf","ovule")
#Sub = c("At","Dt","contig")
## get numbers of each tissue
getsubgenome = function(circbase,Tissue){
  circ.At = filter(circbase,str_detect(circbase$chr,"^A"))
  circ.Dt = filter(circbase,str_detect(circbase$chr, "^D"))
  circ.contig = filter(circbase,str_detect(circbase$chr,"^C"))
  x = data.frame(Tissue = Tissue,
                 Subgenome = c("At","Dt","Contig"),
                 Num = 0)
  x[1,3] = dim(circ.At)[1]
  x[2,3] = dim(circ.Dt)[1]
  x[3,3] = dim(circ.contig)[1]
  return(x)
}
fiberNum = getsubgenome(tmp.F,"fiber")
leafNum = getsubgenome(tmp.L,"leaf")
ovuleNum = getsubgenome(tmp.O,"ovule")
## merge table 
data.plot = rbind(fiberNum,leafNum,ovuleNum)
# data.plot

data.plot$Subgenome = factor(data.plot$Subgenome,levels = c("At","Dt","Contig"))
data.plot$Tissue = factor(data.plot$Tissue,levels = c("fiber","leaf","ovule"))
data.plot$postion = c(228,114,14,178,86,3,295,114,10)
P  = ggplot(data.plot,aes(x = Tissue,y = Num, fill = Subgenome))+
  geom_bar(stat ="identity",width = 0,position =position_stack(reverse = F))+
  scale_fill_brewer(palette = 'Accent')+
  labs(y = "CircRNA number")+
  theme_bw()+
  geom_text(aes(label = data.plot$Num,y = data.plot$postion),vjust = -0,size = 3,hjust = 0)+
  theme(legend.position = "top",axis.title.x = element_blank(),
        panel.grid = element_line(linetype = "dashed",colour = "grey",size = 0),
        panel.border = element_rect(size = 1),
        strip.background = element_rect(size = 1,fill = "#FAFAD2"),
        strip.text = element_text(face = "bold.italic"))
ggsave(filename = "~/03.Project/03.circRNA/CircRNAVersion3/04.result/01.Indentification/CircNumofSubgenomeALL.pdf",plot = P,
       width = 7,height = 5,dpi = 300)
ggsave(filename = "~/03.Project/03.circRNA/CircRNAVersion3/04.result/01.Indentification/CircNumofSubgenomeALL.svg",plot = P,
       width = 7,height = 5,dpi = 300)
ggsave(filename = "~/03.Project/03.circRNA/CircRNAVersion3/04.result/01.Indentification/CircNumofSubgenomeALL.tiff",plot = P,
       width = 7,height = 5,dpi = 300)
View(circNum)
## export circIDs

out = function(path,name,data){
  write.table(x = data.frame(circID = data[,1]),
              file = paste(path,name,".xls",sep = ""),
              sep = "\t",
              row.names = F,
              quote = F,col.names = F)
}
out(data = fiber1,path = "~/03.Project/03.circRNA/CircRNAVersion3/03.process/03.Subgenome/",name = "fiber")
out(data = leaf1,path = "~/03.Project/03.circRNA/CircRNAVersion3/03.process/03.Subgenome/",name = "leaf")
out(data = ovule1,path = "~/03.Project/03.circRNA/CircRNAVersion3/03.process/03.Subgenome/",name = "ovule")

## merge
x.fiber = data.frame(circID = fiber1$circID,
                     tissue = "fiber")
x.ovule = data.frame(circID = ovule1$circID,
                     tissue = "ovule")
x.leaf = data.frame(circID = leaf1$circID,
                    tissue = "leaf")
x = rbind(x.fiber,x.ovule,x.leaf)
head(x)
x$tissue = factor(x$tissue,levels = c("fiber","leaf","ovule"))
save(x,file = "~/03.Project/03.circRNA/CircRNAVersion3/02.data/03.workingdir/BtJ-upset.circID.Rdata")
