#######################################################
#     Prj: Heterosis-circ
#     Assignment: distribution of circRNA on gh chr
#     Author: Shawn Wang
#     Date: Sep 09, 2020
#######################################################
setwd("~/03.Project/03.circRNA/CircRNAVersion3/03.process/03.Subgenome/")
table = read.table("~/03.Project/03.circRNA/CircRNAVersion3/02.data/03.workingdir/chrdistr.table",
                   header = T,
                   sep = "\t")
table = arrange(table,chr)
head(table)
at = filter(table,str_detect(string = table$chr,pattern = "A"))
dt = filter(table,str_detect(string = table$chr,pattern = "D"))
contig = filter(table,str_detect(string = table$chr,pattern = "C"))
table = rbind(at,dt,contig)
table$Subgenome = factor(c(rep("At",times = 13),rep("Dt",times = 13),"Contig"),levels = c("At","Dt","Contig"))
table$chr = factor(table$chr,levels = table$chr)
p =ggplot(data = table,mapping = aes(x = chr,y = number,fill = Subgenome)) +
  geom_bar(stat = "identity") + 
  scale_fill_brewer(palette = 'Accent')+
  labs(y = "CircRNA number") +
  theme_bw()+
  geom_text(aes(label = table$number),size = 3, vjust = -0.2,hjust = 0.5)+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5),legend.position = "top",axis.title.x = element_blank(),
        panel.grid = element_line(linetype = "dashed",colour = "grey",size = 0.1),
        panel.border = element_rect(size = 1),
        strip.background = element_rect(size = 1,fill = "#FAFAD2"),
        strip.text = element_text(face = "bold.italic"))
ggsave(filename = "~/03.Project/D",plot = p,
       width = 7,height = 5,dpi = 300)

# p 
