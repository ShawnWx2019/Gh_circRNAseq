###############################################################
#     Prj: Heterosis-circRNA
#     Assignment: circRNA and hosts Structure visualization
#     Date: Nov 10, 2020
#     Author: Shawn Wang
###############################################################
library(ggplot2)
library(dplyr)
options(stringsAsFactors = F)
setwd("~/03.project/03.circRNA/CircRNAVersion3/03.process/06.SPENAE/structure/")
drawGeneStructure = function(gff,size,lcolor){
  names(gff) = c("chr","source","type","start","end","score","strand","phase","attributes")
  chr = unique(gff$chr)
  strand = unique(gff$strand)
  a = filter(gff,type == "gene")$attributes
  desc = strsplit(a,split = ";")%>%unlist(.)
  tpm1 = filter(gff,type == "exon")
  tpm3 = filter(gff,type == "five_prime_UTR"| type == "three_prime_UTR")
  if (strand == "-"){
    tpm2 = data.frame(type = paste(tpm1$type,seq(1:nrow(tpm1)),sep = ""),
                      end = tpm1$start,
                      start = tpm1$end,
                      y = 1)
    xend = tpm2$end[1]-100
  } else {
    tpm2 = data.frame(type = paste(tpm1$type,seq(1:nrow(tpm1)),sep = ""),
                      start = tpm1$start,
                      end = tpm1$end,
                      y = 1)
    xend = tpm2$end[1]+100
  }
  p = ggplot(tpm2,aes(x = start,y = y))+
    geom_segment(aes(x = tpm2$start[1],y = 0, xend = tpm2$end[nrow(tpm2)], yend = 0),
                 color = lcolor)+
    geom_segment(data = tpm2,
                 mapping = aes(x = start,xend = end,y = 0, yend = 0,color = type),
                 size = size)+
    geom_segment(aes(x = tpm2$end[1],y = 0, xend = xend, yend = 0),
                 arrow = arrow(length = unit(0,"cm")),
                 color = "red")+
    xlab(chr)+
    theme(axis.ticks.y.left = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          axis.line.x.top = element_blank(),
          panel.grid = element_blank(),
          axis.text.y = element_blank(),
          plot.background = element_blank(),
          panel.background = element_blank(),
          axis.title.y = element_blank(),
          axis.line = element_line(colour = "black",size = 1),
          legend.position = "none")
  
  if(nrow(tpm3) != 0){
    p = p +
      geom_segment(data = tpm3,
                   mapping = aes(x = start,xend = end,y = 0, yend = 0),
                   size = size,
                   color = "grey")
  } else {
    p = p
  }
  return(p)
}
circplot = function(gff,order,color,size){
  names(gff) = c("chr","source","type","start","end","score","strand","phase","attributes")
  chr = unique(gff$chr)
  strand = unique(gff$strand)
  tpm1 = filter(gff,type == "exon")
  if (strand == "-") {
    a = (nrow(tpm1)+1)-order
    x = tpm1[a,c(4,5)]
    len <- 0
    y <- 0
    for (i in 1:nrow(x)) {
      len[i] = x[i,2]-x[i,1]
      y[i+1] = y[i] + len[i]
    }
    x = data.frame(start = y[-length(y)],
                   end = y[-1])
    p = ggplot()+
      geom_segment(data = x,
                   mapping = aes(x = start,xend = end,y = 0, yend = 0),
                   color = color,
                   size = size)+ 
      coord_polar()+
      theme_void()
  }
}
drawpic = function(path,name,size,lcolor){
  gff = read.delim(path,
                   header = F,
                   sep = "\t")
  p = drawGeneStructure(gff = gff,size = 4,lcolor = "black")
  ggsave(plot = p,filename = paste(name,".pdf",sep = ""), width = 8,height = 1,dpi = 300)
  assign("gff",value = gff, envir = globalenv())
}
dev.off
## for circ0968
drawpic(path = "968-Gh_A08G039000.txt",name ="968-Gh_A08G039000",size = 4,lcolor = "black" )
p2 = circplot(gff = gff,order = c(3,4),color = c("#e76bf3","#9590ff"),size = 4)
ggsave(plot = p2,filename = "circ_0000968.pdf", width = 2,height = 2,dpi = 300)
## for circ0657
drawpic(path = "657-Gh_A06G235600.txt",name ="657-Gh_A06G235600",size = 4,lcolor = "black" )
p2 = circplot(gff = gff,order = c(5,6),color = c("#53b400","#c49a00"),size = 4)
ggsave(plot = p2,filename = "circ_0000657.pdf",width = 2,height = 2,dpi = 300)
## for circ0675
drawpic(path = "675-Gh_A06G075200.txt",name = "675-Gh_A06G075200",size = 4, lcolor = "black")
p2 = circplot(gff = gff,order = c(5,6),color = c("#c09bff","#00a6ff"),size = 4)
ggsave(plot = p2,filename = "circ_0000675.pdf",width = 2,height = 2,dpi = 300)
## for circ0099
drawpic(path = "99-Gh_A01G160900.txt",name = "99-Gh_A01G160900",size = 4, lcolor = "black")
p2 = circplot(gff = gff,order = c(4,6,7),color = c("#7cae00","#cd9600","#e68613"),size = 4)
ggsave(plot = p2,filename = "circ_0000099.pdf",width = 2,height = 2,dpi = 300)
## for circ0051-- Intronic circRNA
drawpic(path = "51-Gh_A01G114000.txt",name = "51-Gh_A01G114000",size = 4, lcolor = "black")
gff = read.delim("51-Gh_A01G114000.txt",
                 header = F,
                 sep = "\t")
p = ggplot()+
  geom_segment(mapping = aes(x = 21405545,xend = 21405864,y = 0, yend = 0),
               color = "black",
               size = 4)+
  coord_polar()+
  theme_void()
ggsave(plot = p,filename = "circ_0000051.pdf",width = 2,height = 2,dpi = 300)
## for circ2575
drawpic(path = "2575-Gh_D06G073500.txt",name = "2575-Gh_D06G073500",size = 4, lcolor = "black")
p2 = circplot(gff = gff,order = c(4,5),color = c("#9590ff","#00b0f6"),size = 4)
ggsave(plot = p2,filename = "circ_00002575.pdf",width = 2,height = 2,dpi = 300)
## for circ1675
circ1675 = read.delim("~/03.project/03.circRNA/CircRNAVersion3/03.process/06.SPENAE/circ_1675.txt",
                      header = T,
                      sep = "\t",
                      stringsAsFactors = F)

circbody = filter(circ1675,feature == "circ")
genebody = filter(circ1675,feature != "circ")
p = ggplot(circ1675,aes(x = start,y = 1))+
  geom_segment(aes(x = circ1675$start[1],y = 0, xend = circ1675$end[nrow(circ1675)], yend = 0),
               size = 1, color = "black")+
  geom_segment(data = circbody,
               mapping = aes(x = start,xend = end,y = 0, yend = 0),
               size = 4,
               color = c("purple","green","gold"))+
  geom_segment(data = genebody,
               mapping = aes(x = start,xend = end,y = 0, yend = 0),
               size = 2,
               color = "salmon")+
  xlab("A13")+
  theme(axis.ticks.y.left = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x.top = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(colour = "black",size = 1),
        legend.position = "none")
ggsave(plot = p,filename = "structure-1675.pdf",width = 8,height = 1,dpi = 300)

x = circbody[,c(2,3)]
len <- 0
y <- 0
for (i in 1:nrow(x)) {
  len[i] = x[i,2]-x[i,1]
  y[i+1] = y[i] + len[i]
}
x = data.frame(start = y[-length(y)],
               end = y[-1])
p = ggplot()+
  geom_segment(data =x,
               mapping = aes(x = start,xend = end,y = 0, yend = 0),
               color = c("purple","green","gold"),
               size = 4)+ 
  coord_polar()+
  theme_void()
ggsave(plot = p,filename = "circ00001675.pdf",width = 2,height = 2,dpi = 300)
