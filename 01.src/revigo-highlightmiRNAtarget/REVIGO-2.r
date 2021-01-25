setwd("/Volumes/Samsung_T5/03.circRNA/CircRNAVersion3/03.process/07.Sponge/")

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
library(ggrepel)

# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0000003","reproduction", 0.769,-1.505, 0.249, 4.994,-1.6399,1.000,0.000),
c("GO:0001666","response to hypoxia", 0.049,-2.663,-4.660, 3.802,-2.8616,0.949,0.000),
c("GO:0006184","(obsolete) GTP catabolic process", 0.143,-5.371,-1.575, 4.254,-2.3695,0.993,0.000),
c("GO:0006200","(obsolete) ATP catabolic process", 0.143,-4.310,-2.058, 4.254,-3.6925,0.993,0.000),
c("GO:0007155","cell adhesion", 0.544,-3.064, 1.106, 4.844,-5.1002,0.993,0.000),
c("GO:0022610","biological adhesion", 0.550,-2.280, 4.947, 4.849,-5.1002,0.993,0.000),
c("GO:0032501","multicellular organismal process", 2.373,-4.247, 1.165, 5.483,-1.6055,0.993,0.000),
c("GO:0032502","developmental process", 2.812,-3.677, 1.157, 5.557,-1.4195,0.993,0.000),
c("GO:0040007","growth", 0.317, 0.597, 0.783, 4.609,-4.9208,0.993,0.000),
c("GO:0060151","peroxisome localization", 0.001,-5.107, 3.983, 1.903,-2.8428,0.946,0.000),
c("GO:1901990","regulation of mitotic cell cycle phase transition", 0.158, 1.656,-6.954, 4.308,-6.9101,0.726,0.000),
c("GO:1905393","plant organ formation", 0.009, 3.302, 6.165, 3.054,-6.3556,0.827,0.000),
c("GO:0006513","protein monoubiquitination", 0.016,-0.058, 2.635, 3.318,-4.9586,0.912,0.021),
c("GO:0009812","flavonoid metabolic process", 0.018, 1.515, 2.424, 3.355,-1.6006,0.977,0.026),
c("GO:0071554","cell wall organization or biogenesis", 0.950, 0.976, 3.887, 5.086,-3.2628,0.961,0.029),
c("GO:0050665","hydrogen peroxide biosynthetic process", 0.003,-2.482, 0.661, 2.526,-3.4078,0.946,0.031),
c("GO:0051186","cofactor metabolic process", 3.985, 0.686, 0.082, 5.709,-2.6594,0.942,0.050),
c("GO:0016125","sterol metabolic process", 0.106, 5.883, 1.093, 4.135,-3.2782,0.861,0.061),
c("GO:0072593","reactive oxygen species metabolic process", 0.282,-0.614, 1.307, 4.558,-1.3181,0.953,0.065),
c("GO:0009813","flavonoid biosynthetic process", 0.016,-1.854,-1.387, 3.317,-1.8072,0.965,0.082),
c("GO:0017000","antibiotic biosynthetic process", 0.028,-0.064,-2.120, 3.560,-1.3362,0.947,0.093),
c("GO:0042732","D-xylose metabolic process", 0.014, 5.537, 1.722, 3.253,-2.5125,0.897,0.100),
c("GO:0046274","lignin catabolic process", 0.010, 3.525,-2.354, 3.116,-2.1702,0.859,0.119),
c("GO:0016049","cell growth", 0.153, 7.136,-4.077, 4.294,-4.9318,0.860,0.145),
c("GO:0006928","movement of cell or subcellular component", 0.973, 5.823,-2.782, 5.097,-3.8729,0.860,0.169),
c("GO:0007018","microtubule-based movement", 0.287, 5.370,-3.393, 4.567,-2.7529,0.851,0.178),
c("GO:0030029","actin filament-based process", 0.398, 5.612,-2.810, 4.708,-3.3747,0.870,0.184),
c("GO:0009914","hormone transport", 0.074,-6.047,-3.133, 3.980,-2.5441,0.788,0.185),
c("GO:0007059","chromosome segregation", 0.476, 6.145,-3.812, 4.786,-3.7077,0.868,0.187),
c("GO:0007017","microtubule-based process", 0.658, 6.276,-3.335, 4.927,-3.0177,0.865,0.194),
c("GO:0006275","regulation of DNA replication", 0.116,-0.712,-7.572, 4.172,-2.0635,0.810,0.205),
c("GO:0044089","positive regulation of cellular component biogenesis", 0.193,-0.288,-7.283, 4.394,-4.4841,0.675,0.208),
c("GO:0032970","regulation of actin filament-based process", 0.179, 1.837,-7.117, 4.362,-3.9547,0.747,0.212),
c("GO:0046653","tetrahydrofolate metabolic process", 0.251, 6.792,-2.687, 4.508,-2.4529,0.819,0.216),
c("GO:0098661","inorganic anion transmembrane transport", 0.470,-6.659, 2.548, 4.780,-2.5244,0.955,0.228),
c("GO:0016197","endosomal transport", 0.131,-6.734, 2.104, 4.225,-2.2271,0.958,0.238),
c("GO:0007034","vacuolar transport", 0.133,-6.185, 2.556, 4.231,-1.4844,0.958,0.239),
c("GO:0042558","pteridine-containing compound metabolic process", 0.497, 3.436,-2.671, 4.804,-1.4082,0.920,0.250),
c("GO:0048193","Golgi vesicle transport", 0.297,-5.873, 2.687, 4.581,-1.6353,0.956,0.256),
c("GO:0048518","positive regulation of biological process", 1.744,-2.096,-7.873, 5.350,-2.3855,0.848,0.257),
c("GO:0006730","one-carbon metabolic process", 0.328, 7.511,-2.013, 4.625,-2.0672,0.845,0.274),
c("GO:0043085","positive regulation of catalytic activity", 0.818,-2.253,-7.706, 5.021,-2.2701,0.819,0.285),
c("GO:0016192","vesicle-mediated transport", 1.085,-7.021, 2.476, 5.144,-1.3718,0.952,0.289),
c("GO:0046777","protein autophosphorylation", 0.077,-0.204, 2.534, 3.993,-1.3101,0.888,0.297),
c("GO:0009966","regulation of signal transduction", 0.857,-2.079,-6.666, 5.041,-2.3516,0.754,0.300),
c("GO:0009817","defense response to fungus, incompatible interaction", 0.002,-3.306,-4.719, 2.501,-2.4669,0.963,0.323),
c("GO:0065008","regulation of biological quality", 3.395,-1.850,-7.561, 5.639,-1.8469,0.845,0.339),
c("GO:0016122","xanthophyll metabolic process", 0.002, 6.688, 1.109, 2.365,-3.0691,0.877,0.363),
c("GO:0006544","glycine metabolic process", 0.221, 7.200,-1.912, 4.452,-2.2420,0.836,0.366),
c("GO:0006611","protein export from nucleus", 0.020,-6.263, 3.707, 3.399,-2.0175,0.917,0.371),
c("GO:0042133","neurotransmitter metabolic process", 0.006,-3.718,-6.935, 2.903,-1.5023,0.863,0.387),
c("GO:0031408","oxylipin biosynthetic process", 0.007, 7.208,-0.085, 2.941,-1.7987,0.844,0.393),
c("GO:0009734","auxin-activated signaling pathway", 0.035,-2.113,-6.453, 3.649,-2.0317,0.804,0.398),
c("GO:0007010","cytoskeleton organization", 0.786, 2.354,-1.639, 5.004,-3.1599,0.803,0.405),
c("GO:0046688","response to copper ion", 0.020,-2.893,-4.987, 3.403,-1.8652,0.958,0.409),
c("GO:0005985","sucrose metabolic process", 0.017, 4.780,-0.984, 3.327,-1.9920,0.873,0.411),
c("GO:0048444","floral organ morphogenesis", 0.003, 2.786, 6.555, 2.582,-5.8508,0.743,0.416),
c("GO:0009606","tropism", 0.006,-4.513,-4.717, 2.884,-1.4459,0.963,0.416),
c("GO:0010090","trichome morphogenesis", 0.003, 4.761, 4.800, 2.617,-4.0044,0.676,0.417),
c("GO:0051028","mRNA transport", 0.075,-6.400, 3.284, 3.986,-1.8632,0.924,0.431),
c("GO:0009880","embryonic pattern specification", 0.025, 2.693, 6.320, 3.513,-2.9115,0.769,0.438),
c("GO:0001505","regulation of neurotransmitter levels", 0.055,-3.573,-7.151, 3.848,-1.5023,0.866,0.444),
c("GO:0090627","plant epidermal cell differentiation", 0.007, 4.373, 5.127, 2.958,-1.4650,0.764,0.445),
c("GO:0021700","developmental maturation", 0.074, 3.624, 6.351, 3.977,-2.0725,0.783,0.454),
c("GO:0098727","maintenance of cell number", 0.036, 3.067, 6.107, 3.660,-2.0154,0.826,0.458),
c("GO:0097435","supramolecular fiber organization", 0.345, 5.965,-4.770, 4.646,-2.9371,0.767,0.463),
c("GO:0031503","protein complex localization", 0.062,-6.113, 3.703, 3.899,-1.7922,0.942,0.465),
c("GO:0044706","multi-multicellular organism process", 0.067, 0.848, 7.101, 3.933,-1.4779,0.895,0.469),
c("GO:0019827","stem cell population maintenance", 0.035, 3.237, 6.587, 3.651,-2.0154,0.773,0.484),
c("GO:0019439","aromatic compound catabolic process", 1.164, 2.131,-4.345, 5.174,-1.5927,0.919,0.485),
c("GO:0006464","cellular protein modification process", 7.726,-0.575,-1.285, 5.996,-1.3142,0.876,0.485),
c("GO:0009555","pollen development", 0.019, 3.218, 6.584, 3.382,-1.3319,0.779,0.486),
c("GO:0008202","steroid metabolic process", 0.161, 5.375, 1.024, 4.315,-2.7492,0.878,0.486),
c("GO:0006403","RNA localization", 0.118,-6.024, 3.372, 4.179,-1.5969,0.940,0.487),
c("GO:0090066","regulation of anatomical structure size", 0.216,-3.252,-7.386, 4.444,-1.8731,0.853,0.490),
c("GO:0031407","oxylipin metabolic process", 0.007, 7.211,-0.015, 2.961,-1.3711,0.853,0.490),
c("GO:0071166","ribonucleoprotein complex localization", 0.097,-6.702, 3.551, 4.094,-1.7922,0.934,0.492));

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

ggsave("./sim.BP.pdf",plot = p1,width = 10,height = 10,dpi = 300);

