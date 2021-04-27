


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
library(ggrepel);

# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0003824","catalytic activity",6527, 363, 464, 667,-2975,093,0.000),
c("GO:0004871","signal transducer activity", 278,-082, 1.079, 592,-1320,080,0.000),
c("GO:0004888","transmembrane signaling receptor activity", 062, 126,-461, 532,-1220,036,0.000),
c("GO:0005488","binding",5556, 246, 365, 694,-2860,091,0.000),
c("GO:0008144","drug binding", 078, 248, 7.017, 400,-8464,044,0.000),
c("GO:0015562","efflux transmembrane transporter activity", 0.005,-134,-055, 212,-2963,014,0.000),
c("GO:0016832","aldehyde-lyase activity", 001,-089,-194, 452,-4591,052,0.000),
c("GO:0030695","GTPase regulator activity", 006, 487, 691, 463,-1822,063,0.000),
c("GO:0052716","hydroquinone:oxygen oxidoreductase activity", 0.010,-769, 007, 342,-2041,023,0.019),
c("GO:0004683","calmodulin-dependent protein kinase activity", 0.023, 714, 100, 307,-3327,054,0.020),
c("GO:0003756","protein disulfide isomerase activity", 0.025,-517, 023, 352,-1248,034,0.020),
c("GO:0003774","motor activity", 099,-291, 640, 450,-4556,099,0.024),
c("GO:0016787","hydrolase activity",2294, 2.052,-202, 697,-2428,037,0.041),
c("GO:0008092","cytoskeletal protein binding", 008,-242, 146, 499,-4467,098,0.046),
c("GO:0030246","carbohydrate binding", 023, 440,-602, 5.007,-1257,038,0.053),
c("GO:0005525","GTP binding", 183,-481,-300, 500,-1713,007,0.058),
c("GO:0044877","macromolecular complex binding", 040, 310,-680, 5.018,-1320,037,0.058),
c("GO:0097367","carbohydrate derivative binding",1752,-126,-518, 685,-7601,009,0.090),
c("GO:0016740","transferase activity",21.036, 054, 488, 672,-2.0387,038,0.093),
c("GO:0016210","naringenin-chalcone synthase activity", 0.001, 485,-0.082, 274,-3478,009,026),
c("GO:0004372","glycine hydroxymethyltransferase activity", 0.050, 7.037,-299, 346,-2594,092,054),
c("GO:0004568","chitinase activity", 0.021,-235, 456, 378,-3702,082,067),
c("GO:0043167","ion binding",3392,-037,-7.018, 673,-6126,097,082),
c("GO:0004842","ubiquitin-protein transferase activity", 052, 527, 083, 495,-1422,085,086),
c("GO:0004620","phospholipase activity", 030,-363, 568, 463,-2476,046,092),
c("GO:1901265","nucleoside phosphate binding",2085,-355,-544, 654,-5190,064,093),
c("GO:0032440","2-alkenal reductase [NAD(P)] activity", 0.012,-709, 047, 324,-1929,047,096),
c("GO:0036094","small molecule binding",2137,-010,-663, 678,-6556,005,097),
c("GO:0019787","ubiquitin-like protein transferase activity", 078, 653, 042, 426,-1706,084,019),
c("GO:0008238","exopeptidase activity", 065,-284, 519, 5.086,-3056,008,027),
c("GO:0017171","serine hydrolase activity", 138,-156, 665, 541,-1709,042,055),
c("GO:0043168","anion binding",2042,-319,-566, 670,-8319,092,060),
c("GO:0097159","organic cyclic compound binding",4137,-137,-600, 663,-3536,093,062),
c("GO:1901363","heterocyclic compound binding",4115,-244,-612, 663,-3289,093,092),
c("GO:0016772","transferase activity, transferring phosphorus-containing groups", 804, 744,-016, 6.078,-3101,050,006),
c("GO:0070566","adenylyltransferase activity", 057, 745, 092, 459,-3546,052,014),
c("GO:0016817","hydrolase activity, acting on acid anhydrides", 723,-2.063, 649, 6.007,-5899,013,032),
c("GO:0019001","guanyl nucleotide binding", 146,-5.098,-3.007, 515,-1613,030,082),
c("GO:0004712","protein serine/threonine/tyrosine kinase activity", 0.028, 677,-105, 396,-2843,057,095),
c("GO:0003777","microtubule motor activity", 018,-279, 664, 487,-2874,009,031),
c("GO:0035639","purine ribonucleoside triphosphate binding",1515,-456,-425, 648,-9518,095,038),
c("GO:0016298","lipase activity", 087,-4.059, 524, 421,-2641,049,053),
c("GO:0016779","nucleotidyltransferase activity", 154, 7.035, 0.078, 539,-2523,027,056),
c("GO:0016905","myosin heavy chain kinase activity", 0.001, 634,-198, 276,-2523,073,065),
c("GO:0016407","acetyltransferase activity", 131, 688,-025, 539,-1199,067,073),
c("GO:0016857","racemase and epimerase activity, acting on carbohydrates and derivatives", 038,-535, 077, 425,-1094,033,091),
c("GO:0016742","hydroxymethyl-, formyl- and related transferase activity", 032, 661, 143, 414,-1604,082,095));

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
ggsave("./sim.MF.pdf",plot = p1,width = 10,height = 10,dpi = 300);
