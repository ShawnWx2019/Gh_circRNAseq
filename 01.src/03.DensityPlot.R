#####################################
#     Prj: Heterosis-circ
#     Assignment: Re-draw density plot
#     Author: Shawn Wang
#     Date: Sep 09, 2020
####################################
options(stringsAsFactors = F)
library(stringr)
library(dplyr)
library(RColorBrewer)
library(cowplot)
library(reshape2)
## 01. table
setwd("~/03.Project/03.circRNA/CircRNAVersion3/03.process/01.Indentification/")
raw = read.delim("~/03.Project/03.circRNA/CircRNAVersion3/02.data/02.CIRIresult/seq_infor/novel_circRNAs.txt",
                 header = T,
                 sep = "\t")
colnames(raw)
clean = raw[,c(1,9)]
head(clean)
exon = data.frame(filter(clean,str_detect(string = clean$feature,pattern = "exon")),
                  Tag = "exon")
intron = data.frame(filter(clean,str_detect(string = clean$feature,pattern = "intron")),
                  Tag = "intron")
intergenic = data.frame(filter(clean,str_detect(string = clean$feature,pattern = "intergenic")),
                        Tag = "intergenic")
categrory = rbind(exon[,-2],intron[,-2],intergenic[,-2])
save(categrory,file = "~/03.Project/03.circRNA/CircRNAVersion3/02.data/03.workingdir/03.CircRNAType.Rdata")
##=======01.Density plot of circRNA Length based on different category of circRNAs====
length = raw[,c("circRNA_id","full_length")]
head(length)
data.len = inner_join(categrory,length,by = "circRNA_id")
head(data.len)
names(data.len) = c("circRNA_id","Type","full_length")
fontsize = 12
## draw main plot
p = ggplot(data = data.len,mapping = aes(x = full_length, fill = Type)) + 
  geom_density(alpha = 0)+
  scale_fill_manual(values = c("#FF69B4","#00FF7F","#FFD700"))+
  geom_vline(aes(xintercept=3000), color = "blue",linetype="dashed")+
  geom_segment(aes(x = 3000,y = 5e-04, xend = 8000, yend = 5e-04),
               arrow = arrow(length = unit(0,"inches")),
               color = "blue",
               linetype = "dashed"
               )+
  theme_linedraw() +
  theme(panel.grid.major=element_line(colour = "grey",linetype = 2, size = 05),
        panel.grid.minor=element_line(colour = "grey",linetype = 2, size = 05))+
  labs(x = "bp",y = "Density",fontsize = fontsize)+
  xlim(0,20000)+
  theme(legend.position = "top")
# draw inset plot
p1 = ggplot(data = data.len,mapping = aes(x = full_length, fill = Type)) + 
  geom_density(alpha = 0)+
  scale_fill_manual(values = c("#FF69B4","#00FF7F","#FFD700"))+
  theme_linedraw() +
  theme(panel.grid.major=element_line(colour = "grey",linetype = 2, size = 05),
        panel.grid.minor=element_line(colour = "grey",linetype = 2, size = 05))+
  labs(x = "bp",y = "Density",fontsize = fontsize)+
  xlim(0,3000)+
  theme(legend.position = "none")

## Merge plots
p.inset = ggdraw()+
  draw_plot(p)+
  draw_plot(p1,x = 07,y = ,width = ,height = )
ggsave(filename = "~/03.Project/03.circRNA/CircRNAVersion3/04.result/01.Indentification/LenDensityPlot.pdf",
       plot = p.inset,
       width = 10,
       height = 7,
       dpi = 300)

