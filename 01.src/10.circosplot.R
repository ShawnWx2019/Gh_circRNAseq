###########################################
#       Prj: Heterosis ==> circRNA
#       Assignment: circos plot
#       Author: Shawn Wang
#       Date: Dec 05,2020
###########################################
setwd("/Volumes/Samsung_T5/03.circRNA/CircRNAVersion3/03.process/05.circos/")
#install.packages("circlize")
library(circlize)
library(stringr)
library(dplyr)
library(stringr)
## 1st ==> need chr length to plot the chr track.
## import chrlen
chrlen <- read.table("/Volumes/Samsung_T5/03.circRNA/CircRNAVersion3/03.process/05.circos/Gh.Chr.Len.xls",
                     header = F,
                     sep = "\t",
                     stringsAsFactors = F)
head(chrlen)
## start from 0
chr = data.frame(Chr = chrlen$V1,
                 start = 0,
                 end = chrlen$V2)
## remove contig
chr = filter(chr,str_detect(chr$Chr,"^[A-D][0-9][0-9]"))
chr$Chr = factor(chr$Chr,levels = chr$Chr)
## 2nd ==> gene density plot. need chr bin with gene number
## import gene bin file
demo <- read.table("/Volumes/Samsung_T5/03.circRNA/CircRNAVersion3/03.process/05.circos/Gh.gene_density.txt",
                   header = F,
                   sep = "\t",
                   stringsAsFactors = F)
demo = data.frame(chr = demo$V1,
                  start = demo$V2 - demo[1,2],
                  end = demo$V2,
                  value = demo$V3)
head(demo)
## 3rd ==> circRNA and miRNA postion
## in this cycle, I'll plot miRNA as blue dot and circRNA as green dot in same layer.
# import miRNA postion file 
mi <- read.delim("/Volumes/Samsung_T5/03.circRNA/CircRNAVersion3/03.process/05.circos/04.mi-pos.xls",
                 header = F,
                 sep = "\t",
                 stringsAsFactors = F)
# add values
mi.pos = data.frame(Chr = mi$V2,
                    start = mi$V3,
                    end = mi$V4,
                    value = 10)
## adjust start and end position
for (i in 1:nrow(mi.pos)) {
  start = mi.pos[i,2]
  end = mi.pos[i,3]
  if (start < end) {
    mi.pos[i,2] = start
    mi.pos[i,3] = end
  } else{
    mi.pos[i,2] = end
    mi.pos[i,3] = start
  }
}
## double check the start-end position
mi.pos[,2] < mi.pos[,3]
head(mi.pos)
## import circ.pos
circ = read.delim("/Volumes/Samsung_T5/03.circRNA/CircRNAVersion3/02.data/03.workingdir/BtJ_circRNAs.txt",
                  header = T,
                  sep = "\t")
head(circ)
circ.clean = circ[,c(1,2,3,4,8)]
head(circ.clean)
circ.pos = data.frame(Chr = circ.clean$chr,
                      start = circ.clean$start,
                      end = circ.clean$end,
                      value = 20)
head(circ.pos)
## merge two position file 
bed1 = circ.pos
bed2 = mi.pos
bed_list = list(bed1,bed2)
## circRNA 2 miRNA
circ2miRNAmap = read.table("/Volumes/Samsung_T5/03.circRNA/CircRNAVersion3/02.data/03.workingdir/circ2miRNA.map",
                           header = T,
                           sep = "\t",
                           stringsAsFactors = F)

head(circ2miRNAmap)
names(circ2miRNAmap) = c("miID","circID")
head(circ.clean)
conn.mi2circ = inner_join(circ.clean,circ2miRNAmap,by = "circID")
head(mi)
## change the upper class colnames into low class
names(mi) = c("miID","Chr","start","end")
mi$miID = gsub(pattern = "MIR",replacement = "miR",mi$miID)

conn.mi2circ2 = inner_join(conn.mi2circ,mi,by = "miID")
head(conn.mi2circ2)
## mi2circ
mi2circ = conn.mi2circ2[,c(2:4,7:9)]
head(mi2circ)
## remove contig
a = dplyr::filter(mi2circ,
              str_detect(mi2circ$chr,"^[A-D][0-9][0-9]"),
              str_detect(mi2circ$Chr,"^[A-D][0-9][0-9]"))
## circPos and mi.Pos
circPos = a[,(1:3)]
mi.Pos = a[,(4:6)]

## miRNA id
miid = data.frame(miID = unique(conn.mi2circ$miID))
mi2mRNA = read.table("/Volumes/Samsung_T5/03.circRNA/CircRNAVersion3/02.data/03.workingdir/mi2mRNA.txt",
                     header = T,
                     sep = "\t",
                     stringsAsFactors = F)
## gene.pos
genepos = read.table("/Volumes/Samsung_T5/03.circRNA/CircRNAVersion3/03.process/05.circos/mRNA.pos",
                     header = F,
                     sep = "\t",
                     stringsAsFactors = F)
names(genepos) = c("GID","chr","start","end")
head(mi2mRNA)
mi2mRNAmap = data.frame(miID = mi2mRNA$miRNA,
                        GID = mi2mRNA$target_gene)
mi2mRNAmap$GID = gsub(pattern = "..$",replacement = "",mi2mRNAmap$GID)
btjmi2mRNA = inner_join(miid,mi2mRNAmap,by = "miID")
btjmi2mRNA = unique(btjmi2mRNA)
head(btjmi2mRNA)
length(unique(btjmi2mRNA$GID))
length(unique(circ2miRNAmap$circID))
tmp2 = inner_join(btjmi2mRNA,mi,by = "miID")
head(mi)
head(tmp2)
tmp3 = inner_join(tmp2,genepos,by ="GID" )
mi2mRNAall = filter(tmp3,
                    str_detect(tmp3$Chr,"^[A-D][0-9][0-9]"),
                    str_detect(tmp3$chr,"^[A-D][0-9][0-9]"))
# head(tmp3)
# all.miPos = mi2mRNAall[,(3:5)]
# all.mRPos = mi2mRNAall[,(6:8)]

## only plot expattern links
expPattern = read.table("/Volumes/Samsung_T5/03.circRNA/CircRNAVersion3/02.data/03.workingdir/expressionpatternModel",
                        header = T,
                        sep = "\t",
                        stringsAsFactors = F)
exp.mi2circ = inner_join(expPattern,conn.mi2circ2,by = "circID")
############### export list of expPattern-related circ2miRNA relation ################
## for next steps:
# relaCirc2mi = unique(exp.mi2circ[,c(1,2,7)])
# write.table(relaCirc2mi,file = "/Volumes/Samsung_T5/03.circRNA/CircRNAVersion3/03.process/09.finalpathway/expmi2circ.xls",
#             row.names = F,
#             sep = "\t",
#             quote = F)
# head(exp.mi2circ)
# exp.mi2mRNA = inner_join(relaCirc2mi,mi2mRNAmap,by = "miID")
# exp.mi2mRNA = unique(exp.mi2mRNA)
# head(exp.mi2mRNA)
# write.table(exp.mi2mRNA,file = "/Volumes/Samsung_T5/03.circRNA/CircRNAVersion3/03.process/09.finalpathway/expmi2mRNA.xls",
#             row.names = F,
#             sep = "\t",
#             quote = F)
######################################################################################
tmp1 = filter(exp.mi2circ,
              str_detect(exp.mi2circ$chr,"^[A-D][0-9][0-9]"),
              str_detect(exp.mi2circ$Chr,"^[A-D][0-9][0-9]"))
link2.mi = tmp1[,c(8,9,10)]
link2.circ = tmp1[,c(3:5)]
map = data.frame(circID = tmp1$circID,
                 miID = tmp1$miID)
## for mi2mRNA
exp.mi2mRNA = inner_join(map,mi2mRNAmap,"miID")
exp.mi2mRNA = data.frame(GID = exp.mi2mRNA$GID,
                         miID = exp.mi2mRNA$miID)
exp.mi2mRNA = unique(exp.mi2mRNA)
tmp2 = inner_join(exp.mi2mRNA,mi,"miID") %>% inner_join(.,genepos,"GID")
head(tmp2)
mRNA.highlight = data.frame(GID = unique(tmp2$GID)) 
write.table(x = mRNA.highlight,file = "/Volumes/Samsung_T5/03.circRNA/CircRNAVersion3/03.process/05.circos/highlight.mRNA",
            row.names = F,
            sep = "\t",
            quote = F)
tmp3 = filter(tmp2,
              str_detect(tmp2$chr,"^[A-D][0-9][0-9]"),
              str_detect(tmp2$Chr,"^[A-D][0-9][0-9]"))
tmp3 = unique(tmp3)
link1.mi = tmp3[,c(3,4,5)]
link1.mRNA = tmp3[,c(6,7,8)]
##================Part2 plot==========================
circos.clear() ## cleaning last work
## track
circos.genomicInitialize(chr,major.by = 30000000)
## colored At and Dt subgenome with different color
chrcolor = rep(c("turquoise","#FF6347"),each = 13)
circos.genomicTrackPlotRegion(chr,ylim=c(0,0),track.height=0.03,bg.col=chrcolor,cell.padding=c(0.01, 0, 0.01, 0),track.margin=c(0, 0))
## gene density heatmap
min1 = quantile(demo$value)[1]
quant25 = quantile(demo$value)[2]
mid = quantile(demo$value)[3]
quant75 = quantile(demo$value)[4]
max = quantile(demo$value)[5]
col_fun = colorRamp2(c(min1,quant75, max), c("green","black","red"))
circos.genomicTrack(demo, numeric.column = 4,track.height = 0.06,ylim = c(0,1),
                    panel.fun = function(region,value,...){
                      circos.genomicRect(region, value, col = col_fun(value[[1]]), border = NA)
                    })
circos.genomicTrack(bed_list, track.height = 0.06,
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicPoints(region, value, pch = 16, cex = 0, col = i, ...)
                    })
## links

## background2 mi2mRNA
circos.genomicLink(link1.mi, link1.mRNA, col = "#80b1d3", border = "#80b1d3", lwd = 05)
## link3 all mi2circ
circos.genomicLink(circPos, mi.Pos, col = "#bf5b17", border = "#bf5b17", lwd = 05)
## background1 mi2circRNA
circos.genomicLink(link2.mi, link2.circ, col = "#f0027f", border = "#f0027f", lwd = 1)

a = btjmi2mRNA
a = data.frame(a,
               relation = "mi2mRNA")
b = circ2miRNAmap
b = data.frame(b,
               relation = "mi2circRNA")
names(a) = names(b) =c("source","target","relation")
c = rbind(a,b)
d.mRNA = data.frame(ID = unique(a$target),
                    type = "mRNA",
                    label = "",
                    size = 3)
d.circ = data.frame(ID = unique(b$target),
                    type = "circRNA",
                    label = unique(b$target),
                    size = 30)
d.miRNA = data.frame(ID = unique(b$source),
                     type = "miRNA",
                     label = unique(b$source),
                     size = 35)
d = rbind(d.miRNA,d.circ,d.mRNA)
write.table(d, file = "/Volumes/Samsung_T5/03.circRNA/CircRNAVersion3/03.process/07.Sponge/allnode.txt",
            row.names = F,
            quote = F,
            sep = "\t")
write.table(c, file = "/Volumes/Samsung_T5/03.circRNA/CircRNAVersion3/03.process/07.Sponge/all.edge.txt",
            row.names = F,
            quote = F,
            sep = "\t")
