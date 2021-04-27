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


# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0044430","cytoskeletal part", 196,-418, 020, 506,-5868,029,0.000),
c("GO:0099080","supramolecular complex", 040, 408,-147, 426,-4670,046,0.000),
c("GO:0012505","endomembrane system", 211, 052,-696, 542,-2960,049,0.078),
c("GO:0019005","SCF ubiquitin ligase complex", 0.035, 300, 416, 338,-4366,047,008),
c("GO:0030076","light-harvesting complex", 0.035, 257, 668, 338,-2.0885,050,012),
c("GO:0005768","endosome", 019,-760, 182, 497,-4638,084,056),
c("GO:0031982","vesicle", 162,-607, 233, 527,-3325,030,030),
c("GO:0044422","organelle part", 927,-559, 308, 567,-3.0535,021,074),
c("GO:0005739","mitochondrion", 256,-454, 238, 527,-1011,081,007),
c("GO:0031301","integral component of organelle membrane", 014,-553,-014, 423,-3547,067,008),
c("GO:0005759","mitochondrial matrix", 036,-592, 1.006, 420,-2059,068,038),
c("GO:0009706","chloroplast inner membrane", 0.006,-554,-176, 240,-2110,024,081));

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
ggsave("./sim.CC.pdf",plot = p1,width = 10,height = 10,dpi = 300);
# ggsave("C:/Users/path_to_your_file/revigo-plot.pdf");
