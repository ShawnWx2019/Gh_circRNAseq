#######################################
#     Prj: Heterosis-circ
#     Assignment: WGCNA-Tissue sepcific
#     Author: Shawn Wang
#     Date: Oct 30, 2020
#######################################
## based on my OneStepWGCNA : https://gitee.com/shawnmagic/OneStepWGCNA
## Considering the small number of circRNAs, I modified some parameters to make WGCNA result more reliable and more robust.
setwd("~/03.Project/03.circRNA/CircRNAVersion3/03.process/01.Indentification/01.WGCNA/")
## preparetion
library(DESeq2)
library(ggplot2)
library(dplyr)
library(WGCNA)
library(stringr)
library(ape)
library(ggpubr)
library(pheatmap)
options(stringsAsFactors = F)
library(svglite)
## load data
count = read.delim("~/03.Project/03.circRNA/CircRNAVersion3/02.data/02.CIRIresult/quantity/meanscount.txt",
                   header = T,
                   sep = "\t")
BtJ.count = select(count,contains("Bt"),contains("J"),-contains("SY"),-contains("Z12"))
BtJ.count = BtJ.count[,-c(1,5,9,13,17,21,25,29,33)]
BtJ.count = select(BtJ.count,contains("F"),contains("L"),contains("O"))
colnames(BtJ.count)
BtJ.count = data.frame(circID = count$ID,
                       BtJ.count)
head(BtJ.count)
write.table(BtJ.count,file = "~/03.Project/03.circRNA/CircRNAVersion3/02.data/03.workingdir/BtJ.count.xls",
            row.names = F,
            sep = "\t",
            quote = F)
Samples = colnames(BtJ.count)[-1]
tissue = gsub(pattern = ".*_",replacement = "",Samples)%>%gsub(pattern = ".$",replacement = "",.)
tissue = gsub("L","leaf",tissue)%>%gsub("O","ovule",.)%>%gsub("F","fiber",.)
trait= data.frame(sample = Samples,
                  tissue = tissue)
if (ncol(trait) == 2) {
  x <- trait
  Tcol = as.character(unique(x[,2]))
  b <- list()
  for (i in 1:length(Tcol)) {
    b[[i]] = data.frame(row.names = x[,1],
                        levels = ifelse(x[,2] == Tcol[i],1,0))
  }
  c <- bind_cols(b)
  c <- data.frame(row.names = x$name,
                  c)
  colnames(c) = Tcol
  rownames(c) = trait[,1]
  pheTmp <- c
} else {
  pheTmp = data.frame(row.names = trait[,1],
                      trait[,-1])
}
count = data.frame(row.names = BtJ.count[,1],
                   BtJ.count[,-1])

head(count)
condition <- factor(c(rep("case",16),rep("control",11)),levels = c("case","control"))
colData <- data.frame(row.names = colnames(count), condition)
dds <- DESeqDataSetFromMatrix(count, colData, design = ~ condition)
dds <- DESeq(dds)
vsd <- assay(varianceStabilizingTransformation(dds))
datExpr = data.frame(vsd)
dim(datExpr)
type = "unsigned"
corType = "pearson"
corFnc = ifelse(corType=="pearson", cor, bicor)
maxPOutliers = ifelse(corType=="pearson",1,0.05)
robustY = ifelse(corType=="pearson",T,F)
Title = "CircRNA-Tissue-Specific"
## 1.1 functions
source("~/02.MyScript/OneStepWGCNA/01.Rscript/11.02.WGCNA.SFT.R")
source("~/02.MyScript/OneStepWGCNA/01.Rscript/11.03.WGCNA.module.R")
source("~/02.MyScript/OneStepWGCNA/01.Rscript/11.04.WGCNA.moduleTrait.R")
source("~/02.MyScript/OneStepWGCNA/01.Rscript/11.05.WGCNA.HubGene.R")
## 2.2 sft calculate
GNC = 1
head(datExpr)
str(datExpr)
##==================SFT=========================##
WGCNA.SFT(Title = Title,
          datExpr = datExpr,
          GeneNumCut = GNC)
if (is.na(power)){
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                 ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                        ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                               ifelse(type == "unsigned", 6, 12))
                 )
  )
}
net = blockwiseModules(datExpr, power = power, maxBlockSize = nGenes,
                       TOMType = type, minModuleSize = 25,
                       reassignThreshold = 0, mergeCutHeight = 0.6,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType, 
                       maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                       #saveTOMFileBase = paste(Title,".tom",sep = ""),
                       verbose = 3)
table(net$colors)
assign("net",value = net, envir = globalenv())
####============modular construction=================##
# Convert labels to colors for plotting
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
# Plot the dendrogram and the module colors underneath
pdf(file = paste(Title,"Module.pdf",sep = "."),width = 6,height = 6)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

MEs = net$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
Gene2module <- data.frame(GID = colnames(datExpr),
                          Module = moduleColors)
write.table(Gene2module,file = paste(Title,"Gene2module.xls",sep = "_"),
            row.names = F,
            quote = F,
            sep = "\t")
##===============Module-Trait relationship===================##
phenotype = pheTmp
traitData <- phenotype
dim(traitData)
dim(MEs_col)
dim(traitData)
if (corType=="pearson") {
  modTraitCor = cor(MEs_col, traitData, use = "p")
  modTraitP = corPvalueStudent(modTraitCor, nSamples)
} else {
  modTraitCorP = bicorAndPvalue(MEs_col, traitData, robustY=robustY)
  modTraitCor = modTraitCorP$bicor
  modTraitP   = modTraitCorP$p
}

textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
pdf(file = paste(Title,"Module_trait.pdf",sep = "."))
tiff(filename = paste(Title,"Module_trait.tiff",sep = "."),width = 6,height = 6,units="in", compression="lzw", res=300)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData), 
               yLabels = colnames(MEs_col), 
               cex.lab = 0.7, xLabelsAngle = 0, xLabelsAdj = 0.5,
               ySymbols = substr(colnames(MEs_col),3,1000), colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, 
               cex.text = 0.6, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

save(BtJ.count,pheTmp,trait,file = "~/03.Project/03.circRNA/CircRNAVersion3/01.src/02.TissueSpecificWGCNA.Rdata")

## expression each module
trait$sample
colnames(MEs) = colnames(MEs_col)
dataTrait = trait
rownames(dataTrait) = dataTrait$sample
## heatmap
moduleheatmap = function(datExpr,MEs,which.module){
  ME=MEs[,paste("ME",which.module,sep = "")]
  ## heatmap data
  h <- t(scale(datExpr[,moduleColors == which.module]))
  x = reshape2::melt(data = h,id.var = rownames(h))
  ## barplot data
  b <- data.frame(Sample = factor(rownames(datExpr),levels = rownames(datExpr)) ,
                  Eigengene = ME)
  p1 = ggplot(x, aes(Var2, Var1)) + 
    geom_tile(aes(fill = value),colour = alpha("white",0)) + 
    scale_fill_gradient2(low = "green",mid = "black",high = "red")+
    theme_classic()+
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 90,face = "bold",vjust = 0.5,size = 15),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.x = element_blank())
  ## barplot
  p2 = ggplot(b,mapping  = aes(x = Sample,y = Eigengene)) + geom_bar(stat='identity',fill = which.module)+
    theme_classic()+
    theme(axis.title = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  p3 = ggarrange(p1,p2,ncol = 1,nrow = 2)
  return(p3)
}
p3 = moduleheatmap(datExpr = datExpr,MEs=MEs,which.module = "blue")
p3
p4 = moduleheatmap(datExpr = datExpr,MEs = MEs, which.module = "turquoise")
ggsave(filename = "~/03.Project/03.circRNA/CircRNAVersion3/03.process/01.Indentification/01.WGCNA/CircModule.blue.pdf",plot = p3,width = 12,height = 7,dpi = 300)
ggsave(filename = "~/03.Project/03.circRNA/CircRNAVersion3/03.process/01.Indentification/01.WGCNA/CircModule.blue.tiff",plot = p3,width = 12,height = 7,dpi = 300)
ggsave(filename = "~/03.Project/03.circRNA/CircRNAVersion3/03.process/01.Indentification/01.WGCNA/CircModule.turquoise.pdf",plot = p4,width = 12,height = 7,dpi = 300)
ggsave(filename = "~/03.Project/03.circRNA/CircRNAVersion3/03.process/01.Indentification/01.WGCNA/CircModule.turquoise.tiff",plot = p4,width = 12,height = 7,dpi = 300)


dataTrait$num = rep(c(1,2,3),each = 9)
plotMEpairs(datME = MEs,y = dataTrait$num)

connet = abs(cor(datExpr,use = "p"))^6
Alldegrees1=intramodularConnectivity(connet, moduleColors)
which.module = "turquoise"
Ovule = as.data.frame(traitData$ovule)
names(Ovule) = "ovule"
GS1 = as.numeric(cor(Ovule,datExpr,use = "p"))
G.sig = abs(GS1)
color.level = unique(moduleColors)
n = length(color.level)-1
sizeGrWindow(9,6)
par(mfrow = c(2,as.integer(0.5+n/2)))
par(mar = c(4,5,3,1))

for (i in c(1:n)) {
  whichmodule=color.level[[i]];
  restrict1 = (moduleColors==whichmodule);
  verboseScatterplot(Alldegrees1$kWithin[restrict1],
                     G.sig[restrict1], col=moduleColors[restrict1],
                     main=whichmodule,
                     xlab = "Connectivity", ylab = "Gene Significance", abline = TRUE)
}

##

modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
verbosplot = function(tdf,tname,mm){
  names(tdf) = tname
  geneTraitSignificance = as.data.frame(cor(datExpr, tdf, use = "p"));
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
  names(geneTraitSignificance) = paste("GS.", names(tdf), sep="");
  names(GSPvalue) = paste("p.GS.", names(tdf), sep="")
  module = mm
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = paste("CircRNA significance for",tname,sep = " ") ,
                     #main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
}

# leaf = as.data.frame(traitData$leaf)
# p2 = verbosplot(tdf = leaf,tname = "leaf",mm = "blue")

svg(filename = "~/03.Project/03.circRNA/CircRNAVersion3/04.result/01.Indentification/verbosPlotWGCNA.svg",height = 10,width = 7)
#pdf(file = "~/03.Project/03.circRNA/CircRNAVersion3/04.result/01.Indentification/verbosPlotWGCNA.pdf",height = 10,width = 7)
par(mfrow = c(2,1))
ovule = as.data.frame(traitData$ovule);
verbosplot(tdf = ovule,tname = "ovule",mm = "turquoise")
fiber = as.data.frame(traitData$fiber)
verbosplot(tdf = fiber,tname = "fiber",mm = "blue")
dev.off() 
## enrichment analysis
raw = read.delim("~/03.Project/03.circRNA/CircRNAVersion3/02.data/03.workingdir/circ2host",
                     header = T,
                     sep = "\t")
clean = filter(raw,gene_id != "n/a")
head(clean)
m2c = read.delim("CircRNA-Tissue-Specific_Gene2module.xls",
                 header = T,
                 sep = "\t")
head(m2c)
names(m2c) = c("circID","Module")
x = inner_join(m2c,clean,by = "circID")
head(x)
blue = filter(x,Module == "blue")
turquoise = filter(x, Module == "turquoise")
write.table(data.frame(blue = blue[,3]),file = "blue.xls",
            row.names = F,
            sep = "\t",
            quote = F)

write.table(data.frame(turquoise = turquoise[,3]),file = "turquoise.xls",
            row.names = F,
            sep = "\t",
            quote = F)

bluego = read.table("Enrichment/result/blue.xls.GO.Enrichment.final.xls",
                    header = T,
                    sep = "\t")
turquoisego = read.table("Enrichment/result/turquoise.xls.GO.Enrichment.final.xls",
                         header = T,
                         sep = "\t")
godotplot = function(left,right,lname,rname,name){
  ## top 15 terms
  x = left[order(left$P_value)[1:15],]
  y = right[order(right$P_value)[1:15],]
  ## Gene ratio calculate
  x$GeneRatio = x$HitsGenesCountsInSelectedSet/x$AllGenesCountsInSelectedSet
  y$GeneRatio = y$HitsGenesCountsInSelectedSet/y$AllGenesCountsInSelectedSet
  ## ordered the GO name for plot
  x = x[order(x$GeneRatio),]
  y = y[order(y$GeneRatio,decreasing = T),]
  ## data cleanning
  leftdata = data.frame(GO_Name = factor(x$GO_Name,levels = x$GO_Name),
                        P_value = -log10(x$P_value),
                        GeneRatio = x$GeneRatio,
                        Cultivar = rname)
  rightdata = data.frame(GO_Name = factor(y$GO_Name,levels = y$GO_Name),
                         P_value = -log10(y$P_value),
                         GeneRatio = y$GeneRatio,
                         Cultivar = lname)
  ## merge
  z = rbind(leftdata,rightdata)
  ## plot
  p = ggplot(data = z, mapping = aes(x = Cultivar,y = GO_Name))+
    geom_point(aes(size = GeneRatio,color = P_value))+
    scale_color_gradient(low = "blue", high = "red", guide=guide_colorbar(reverse=TRUE))+
    scale_y_discrete(labels=function(y) stringr::str_wrap(y,width=50)) +
    theme_bw() +
    theme(axis.text.x = element_text(colour = "black",
                                     size = 16, vjust =1 ),
          axis.text.y = element_text(colour = "black",
                                     size = 16, hjust =1 ),
          axis.title = element_text(margin=margin(10, 5, 0, 0),
                                    color = "black",size = 18),
          axis.title.y = element_text(angle=90))+ 
    labs(color=expression(-log[10](Qvalue)),
         size="Gene Ratio",
         x = "",
         y = "")
  #return(p)
  ggsave(p,filename = paste("~/03.Project/03.circRNA/CircRNAVersion3/03.process/01.Indentification/01.WGCNA/",name,".svg",sep = ""),width = 9,height = 10,dpi = 300)
}

p = godotplot(left = bluego,right = turquoisego,lname = "blue",rname = "turquoise",name = "GOdotplot")


