


# A plotting R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800


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
revigo.data <- rbind(c("GO:0003824","catalytic activity",65.827, 3.863, 4.964, 6.967,-2.7975,0.993,0.000),
c("GO:0004871","signal transducer activity", 2.778,-0.282, 1.079, 5.592,-1.5320,0.980,0.000),
c("GO:0004888","transmembrane signaling receptor activity", 0.962, 1.726,-4.561, 5.132,-1.9220,0.936,0.000),
c("GO:0005488","binding",55.656, 2.446, 3.865, 6.894,-2.9860,0.991,0.000),
c("GO:0008144","drug binding", 0.178, 2.948, 7.017, 4.400,-8.2464,0.944,0.000),
c("GO:0015562","efflux transmembrane transporter activity", 0.005,-1.734,-0.355, 2.812,-2.7963,0.914,0.000),
c("GO:0016832","aldehyde-lyase activity", 0.201,-0.889,-1.994, 4.452,-4.5591,0.952,0.000),
c("GO:0030695","GTPase regulator activity", 0.206, 4.287, 6.891, 4.463,-1.5822,0.963,0.000),
c("GO:0052716","hydroquinone:oxygen oxidoreductase activity", 0.010,-7.469, 0.307, 3.142,-2.3041,0.923,0.019),
c("GO:0004683","calmodulin-dependent protein kinase activity", 0.023, 7.214, 1.400, 3.507,-3.8327,0.854,0.020),
c("GO:0003756","protein disulfide isomerase activity", 0.025,-5.617, 0.723, 3.552,-1.7248,0.934,0.020),
c("GO:0003774","motor activity", 0.399,-2.491, 6.740, 4.750,-4.3556,0.799,0.024),
c("GO:0016787","hydrolase activity",22.294, 2.052,-2.702, 6.497,-2.7428,0.937,0.041),
c("GO:0008092","cytoskeletal protein binding", 0.708,-2.442, 1.946, 4.999,-4.5467,0.898,0.046),
c("GO:0030246","carbohydrate binding", 0.723, 4.440,-6.502, 5.007,-1.9257,0.938,0.053),
c("GO:0005525","GTP binding", 1.783,-4.781,-3.500, 5.400,-1.6713,0.807,0.058),
c("GO:0044877","macromolecular complex binding", 0.740, 3.210,-6.580, 5.018,-1.5320,0.937,0.058),
c("GO:0097367","carbohydrate derivative binding",17.252,-1.926,-5.618, 6.385,-7.8601,0.909,0.090),
c("GO:0016740","transferase activity",21.036, 0.854, 4.188, 6.472,-2.0387,0.938,0.093),
c("GO:0016210","naringenin-chalcone synthase activity", 0.001, 4.385,-0.082, 2.274,-3.6478,0.909,0.126),
c("GO:0004372","glycine hydroxymethyltransferase activity", 0.050, 7.037,-2.399, 3.846,-2.7594,0.892,0.154),
c("GO:0004568","chitinase activity", 0.021,-2.235, 4.956, 3.478,-3.5702,0.882,0.167),
c("GO:0043167","ion binding",33.492,-0.637,-7.018, 6.673,-6.6126,0.897,0.182),
c("GO:0004842","ubiquitin-protein transferase activity", 0.352, 5.827, 0.483, 4.695,-1.4422,0.885,0.186),
c("GO:0004620","phospholipase activity", 0.130,-3.963, 5.568, 4.263,-2.4476,0.846,0.192),
c("GO:1901265","nucleoside phosphate binding",20.185,-3.855,-5.644, 6.454,-5.7190,0.864,0.193),
c("GO:0032440","2-alkenal reductase [NAD(P)] activity", 0.012,-7.209, 0.947, 3.224,-1.5929,0.947,0.196),
c("GO:0036094","small molecule binding",21.337,-0.510,-6.463, 6.478,-6.3556,0.905,0.197),
c("GO:0019787","ubiquitin-like protein transferase activity", 0.378, 6.553, 0.542, 4.726,-1.3706,0.884,0.219),
c("GO:0008238","exopeptidase activity", 0.865,-2.184, 5.919, 5.086,-3.4056,0.808,0.227),
c("GO:0017171","serine hydrolase activity", 1.238,-1.256, 6.665, 5.241,-1.3709,0.842,0.255),
c("GO:0043168","anion binding",20.942,-3.119,-5.466, 6.470,-8.1319,0.892,0.260),
c("GO:0097159","organic cyclic compound binding",41.137,-1.637,-6.800, 6.763,-3.6536,0.893,0.262),
c("GO:1901363","heterocyclic compound binding",41.115,-2.344,-6.612, 6.763,-3.6289,0.893,0.292),
c("GO:0016772","transferase activity, transferring phosphorus-containing groups", 8.504, 7.144,-0.516, 6.078,-3.9101,0.850,0.306),
c("GO:0070566","adenylyltransferase activity", 0.257, 7.345, 0.492, 4.559,-3.3546,0.852,0.314),
c("GO:0016817","hydrolase activity, acting on acid anhydrides", 7.223,-2.063, 6.449, 6.007,-5.2899,0.813,0.332),
c("GO:0019001","guanyl nucleotide binding", 1.846,-5.098,-3.007, 5.415,-1.6613,0.830,0.382),
c("GO:0004712","protein serine/threonine/tyrosine kinase activity", 0.028, 6.377,-1.905, 3.596,-2.1843,0.857,0.395),
c("GO:0003777","microtubule motor activity", 0.218,-2.979, 6.364, 4.487,-2.1874,0.809,0.431),
c("GO:0035639","purine ribonucleoside triphosphate binding",15.815,-4.156,-4.725, 6.348,-9.2518,0.795,0.438),
c("GO:0016298","lipase activity", 0.187,-4.059, 5.124, 4.421,-2.1641,0.849,0.453),
c("GO:0016779","nucleotidyltransferase activity", 1.954, 7.035, 0.078, 5.439,-2.9523,0.827,0.456),
c("GO:0016905","myosin heavy chain kinase activity", 0.001, 6.134,-1.598, 2.276,-2.7523,0.873,0.465),
c("GO:0016407","acetyltransferase activity", 1.231, 6.688,-0.625, 5.239,-1.3199,0.867,0.473),
c("GO:0016857","racemase and epimerase activity, acting on carbohydrates and derivatives", 0.238,-5.735, 0.377, 4.525,-1.3094,0.933,0.491),
c("GO:0016742","hydroxymethyl-, formyl- and related transferase activity", 0.232, 6.361, 1.843, 4.514,-1.8604,0.882,0.495));

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
p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = log10_p_value, size = plot_size), alpha = I(0.6) ) + scale_size_area();
p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(one.data$log10_p_value), 0) );
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) + scale_size_area();
p1 <- p1 + scale_size( range=c(5, 20)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.15, ]; 
p1 <- p1 + geom_label_repel( data = ex, 
                             aes(plot_X, plot_Y, label = description), 
                             colour = I(alpha("black", 0.85)), 
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
