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


# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0044430","cytoskeletal part", 1.296,-4.218, 0.820, 5.106,-5.4868,0.329,0.000),
c("GO:0099080","supramolecular complex", 0.540, 4.208,-1.947, 4.726,-4.5670,0.946,0.000),
c("GO:0012505","endomembrane system", 2.811, 0.252,-6.296, 5.442,-2.1960,0.849,0.078),
c("GO:0019005","SCF ubiquitin ligase complex", 0.035, 3.600, 4.516, 3.538,-4.2366,0.747,0.108),
c("GO:0030076","light-harvesting complex", 0.035, 2.157, 6.768, 3.538,-2.0885,0.750,0.212),
c("GO:0005768","endosome", 0.319,-7.260, 1.782, 4.497,-4.5638,0.384,0.256),
c("GO:0031982","vesicle", 1.362,-6.307, 2.833, 5.127,-3.3325,0.530,0.330),
c("GO:0044422","organelle part", 9.427,-5.359, 3.508, 5.967,-3.0535,0.521,0.374),
c("GO:0005739","mitochondrion", 2.156,-4.554, 2.538, 5.327,-1.4011,0.481,0.407),
c("GO:0031301","integral component of organelle membrane", 0.214,-5.753,-0.214, 4.323,-3.9547,0.467,0.408),
c("GO:0005759","mitochondrial matrix", 0.336,-5.692, 1.006, 4.520,-2.1059,0.468,0.438),
c("GO:0009706","chloroplast inner membrane", 0.006,-5.654,-1.376, 2.740,-2.6110,0.524,0.481));

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
ggsave("./sim.CC.pdf",plot = p1,width = 10,height = 10,dpi = 300);
# ggsave("C:/Users/path_to_your_file/revigo-plot.pdf");
