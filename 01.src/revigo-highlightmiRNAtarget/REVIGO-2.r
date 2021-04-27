setwd("/Volumes/Samsung_T5/03.circRNA/CircRNAVersion3/03.process/07.Sponge/")

# A plotting R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10371/journal.pone.0021800


# --------------------------------------------------------------------------
# If you don't have the ggplot2 package installed, uncomment the following line:
# install.packages( "ggplot2" );
library( ggplot2 );
# --------------------------------------------------------------------------
# If you don't have the scales package installed, uncomment the following line:
# install.packages( "scales" );
library( scales );
library(ggrepel)

# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0000003","reproduction", 069,-105, 049, 494,-1399,1.000,0.000),
c("GO:0001666","response to hypoxia", 0.049,-263,-460, 302,-2616,049,0.000),
c("GO:0006184","(obsolete) GTP catabolic process", 043,-571,-175, 454,-2695,093,0.000),
c("GO:0006200","(obsolete) ATP catabolic process", 043,-410,-2.058, 454,-3925,093,0.000),
c("GO:0007155","cell adhesion", 044,-3.064, 106, 444,-5002,093,0.000),
c("GO:0022610","biological adhesion", 050,-280, 447, 449,-5002,093,0.000),
c("GO:0032501","multicellular organismal process", 273,-447, 165, 583,-1055,093,0.000),
c("GO:0032502","developmental process", 212,-377, 157, 557,-1195,093,0.000),
c("GO:0040007","growth", 017, 097, 083, 409,-4208,093,0.000),
c("GO:0060151","peroxisome localization", 0.001,-507, 383, 103,-2428,046,0.000),
c("GO:1901990","regulation of mitotic cell cycle phase transition", 058, 156,-654, 408,-6101,026,0.000),
c("GO:1905393","plant organ formation", 0.009, 302, 665, 3.054,-6556,027,0.000),
c("GO:0006513","protein monoubiquitination", 0.016,-0.058, 235, 318,-4586,012,0.021),
c("GO:0009812","flavonoid metabolic process", 0.018, 115, 224, 355,-1006,077,0.026),
c("GO:0071554","cell wall organization or biogenesis", 050, 076, 387, 5.086,-3628,061,0.029),
c("GO:0050665","hydrogen peroxide biosynthetic process", 0.003,-282, 061, 226,-3078,046,0.031),
c("GO:0051186","cofactor metabolic process", 385, 086, 0.082, 509,-2594,042,0.050),
c("GO:0016125","sterol metabolic process", 006, 583, 1.093, 435,-3782,061,0.061),
c("GO:0072593","reactive oxygen species metabolic process", 082,-014, 107, 458,-1181,053,0.065),
c("GO:0009813","flavonoid biosynthetic process", 0.016,-154,-187, 317,-1072,065,0.082),
c("GO:0017000","antibiotic biosynthetic process", 0.028,-0.064,-220, 360,-1362,047,0.093),
c("GO:0042732","D-xylose metabolic process", 0.014, 537, 122, 353,-2125,097,000),
c("GO:0046274","lignin catabolic process", 0.010, 325,-254, 316,-2702,059,019),
c("GO:0016049","cell growth", 053, 736,-4.077, 494,-4318,060,045),
c("GO:0006928","movement of cell or subcellular component", 073, 523,-282, 5.097,-3729,060,069),
c("GO:0007018","microtubule-based movement", 087, 570,-393, 467,-2529,051,078),
c("GO:0030029","actin filament-based process", 098, 512,-210, 408,-3747,070,084),
c("GO:0009914","hormone transport", 0.074,-6.047,-333, 380,-2441,088,085),
c("GO:0007059","chromosome segregation", 076, 645,-312, 486,-3077,068,087),
c("GO:0007017","microtubule-based process", 058, 676,-335, 427,-3.0177,065,094),
c("GO:0006275","regulation of DNA replication", 016,-012,-772, 472,-2.0635,010,005),
c("GO:0044089","positive regulation of cellular component biogenesis", 093,-088,-783, 494,-4841,075,008),
c("GO:0032970","regulation of actin filament-based process", 079, 137,-717, 462,-3547,047,012),
c("GO:0046653","tetrahydrofolate metabolic process", 051, 692,-287, 408,-2529,019,016),
c("GO:0098661","inorganic anion transmembrane transport", 070,-659, 248, 480,-2244,055,028),
c("GO:0016197","endosomal transport", 031,-634, 204, 425,-2271,058,038),
c("GO:0007034","vacuolar transport", 033,-685, 256, 431,-1844,058,039),
c("GO:0042558","pteridine-containing compound metabolic process", 097, 336,-271, 404,-1082,020,050),
c("GO:0048193","Golgi vesicle transport", 097,-573, 287, 481,-1353,056,056),
c("GO:0048518","positive regulation of biological process", 144,-2.096,-773, 550,-2855,048,057),
c("GO:0006730","one-carbon metabolic process", 028, 711,-2.013, 425,-2.0672,045,074),
c("GO:0043085","positive regulation of catalytic activity", 018,-253,-706, 5.021,-2701,019,085),
c("GO:0016192","vesicle-mediated transport", 1.085,-7.021, 276, 544,-1718,052,089),
c("GO:0046777","protein autophosphorylation", 0.077,-004, 234, 393,-1101,088,097),
c("GO:0009966","regulation of signal transduction", 057,-2.079,-666, 5.041,-2516,054,000),
c("GO:0009817","defense response to fungus, incompatible interaction", 0.002,-306,-419, 201,-2669,063,023),
c("GO:0065008","regulation of biological quality", 395,-150,-761, 539,-1469,045,039),
c("GO:0016122","xanthophyll metabolic process", 0.002, 688, 109, 265,-3.0691,077,063),
c("GO:0006544","glycine metabolic process", 021, 700,-112, 452,-2420,036,066),
c("GO:0006611","protein export from nucleus", 0.020,-663, 307, 399,-2.0175,017,071),
c("GO:0042133","neurotransmitter metabolic process", 0.006,-318,-635, 203,-1023,063,087),
c("GO:0031408","oxylipin biosynthetic process", 0.007, 708,-0.085, 241,-1987,044,093),
c("GO:0009734","auxin-activated signaling pathway", 0.035,-213,-653, 349,-2.0317,004,098),
c("GO:0007010","cytoskeleton organization", 086, 254,-139, 5.004,-3599,003,005),
c("GO:0046688","response to copper ion", 0.020,-293,-487, 303,-1652,058,009),
c("GO:0005985","sucrose metabolic process", 0.017, 480,-084, 327,-1920,073,011),
c("GO:0048444","floral organ morphogenesis", 0.003, 286, 655, 282,-5508,043,016),
c("GO:0009606","tropism", 0.006,-413,-417, 284,-1459,063,016),
c("GO:0010090","trichome morphogenesis", 0.003, 461, 400, 217,-4.0044,076,017),
c("GO:0051028","mRNA transport", 0.075,-600, 384, 386,-1632,024,031),
c("GO:0009880","embryonic pattern specification", 0.025, 293, 620, 313,-2115,069,038),
c("GO:0001505","regulation of neurotransmitter levels", 0.055,-373,-751, 348,-1023,066,044),
c("GO:0090627","plant epidermal cell differentiation", 0.007, 473, 527, 258,-1650,064,045),
c("GO:0021700","developmental maturation", 0.074, 324, 651, 377,-2.0725,083,054),
c("GO:0098727","maintenance of cell number", 0.036, 3.067, 607, 360,-2.0154,026,058),
c("GO:0097435","supramolecular fiber organization", 045, 565,-470, 446,-2371,067,063),
c("GO:0031503","protein complex localization", 0.062,-613, 303, 399,-1922,042,065),
c("GO:0044706","multi-multicellular organism process", 0.067, 048, 701, 333,-1779,095,069),
c("GO:0019827","stem cell population maintenance", 0.035, 337, 687, 351,-2.0154,073,084),
c("GO:0019439","aromatic compound catabolic process", 164, 231,-445, 574,-1927,019,085),
c("GO:0006464","cellular protein modification process", 726,-075,-185, 596,-1142,076,085),
c("GO:0009555","pollen development", 0.019, 318, 684, 382,-1319,079,086),
c("GO:0008202","steroid metabolic process", 061, 575, 1.024, 415,-2492,078,086),
c("GO:0006403","RNA localization", 018,-6.024, 372, 479,-1969,040,087),
c("GO:0090066","regulation of anatomical structure size", 016,-352,-786, 444,-1731,053,090),
c("GO:0031407","oxylipin metabolic process", 0.007, 711,-0.015, 261,-1711,053,090),
c("GO:0071166","ribonucleoprotein complex localization", 0.097,-602, 351, 4.094,-1922,034,092));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$plot_size <- as.numeric( as.character(one.data$plot_size) );
one.data$log10_p_value <- as.numeric( as.character(one.data$log10_p_value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );
#head(one.data);


# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

p1 <- ggplot( data = one.data );
p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = log10_p_value, size = plot_size), alpha = I(0) ) + scale_size_area();
p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(one.data$log10_p_value), 0) );
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0) )) + scale_size_area();
p1 <- p1 + scale_size( range=c(5, 20)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 05, ]; 
p1 <- p1 + geom_label_repel( data = ex, 
                             aes(plot_X, plot_Y, label = description), 
                             colour = I(alpha("black", 05)), 
                             arrow = arrow(length = unit(0.01, "npc"), type = "closed", ends = "first"),
                             size = 6,vjust = -3);
p1 <- p1 + labs (y = "semantic space x", x = "semantic space y");
p1 <- p1 + theme(legend.key = element_blank()) ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);

p1 <- p1 + theme(axis.title = element_text(size = 20),
                 axis.text = element_text(size = 18))

# --------------------------------------------------------------------------
# Output the plot to screen

p1;

# Uncomment the line below to also save the plot to a file.
# The file type depends on the extension (default=pdf).

ggsave("./sim.BP.pdf",plot = p1,width = 10,height = 10,dpi = 300);

