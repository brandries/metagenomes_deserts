#Barplots for taxonomic composition
libs <- c("ggplot2", "ggthemes", "extrafont", "plyr", "scales", "wesanderson", "reshape")
lapply(libs, require, character.only = TRUE)

#Set WD
setwd("C:/Users/Andries van der Walt/Dropbox/Masters 2017/Thesis/Data and analysis FC/Figures")

#Import data
#Your data has to be a table with first column being the one of axis, and the second being the other axis.
#the third column contains the values used to draw the bar chart. 
#alternatively, you can input a OTU table style table, and change the format using the reshape::melt function.

compare_genes_func <- read.delim("~/Masters/Data 2017/blast genes/Tsv files/compare_all_genes_norm-L2_filt0.01.txt", sep = "\t", dec = ".", row.names = 1)

#order the figure - pur in the phyla in your own dataset
corr_order <- c("Other", "Cell Growth and Death, Cell Motility", "Infectious Diseases", 
                "Biosynthesis of Other Secondary Metabolites", "Glycan Biosynthesis and Metabolism", 
                "Metabolism of Terpenoids and Polyketides", "Folding", "Sorting and Degradation", "Signal Transduction", 
                "Xenobiotics Biodegradation and Metabolism", "Metabolism of Other Amino Acid", "Lipid Metabolism", 
                "Replication and Repair", "Translation", "Membrane Transport", "Nucleotide Metabolism", 
                "Metabolism of Cofactors and Vitamins", "Energy Metabolism", "Amino Acid Metabolism", 
                "Carbohydrate Metabolism")
corr_nam <- c("Other", "Cell Growth and Death, Cell Motility", "Infectious Diseases", 
              "Biosynthesis of Other Secondary Metabolites", "Glycan Biosynthesis and Metabolism", 
              "Metabolism of Terpenoids and Polyketides", "Folding", "Sorting and Degradation", "Signal Transduction", 
              "Xenobiotics Biodegradation and Metabolism", "Metabolism of Other Amino Acid", "Lipid Metabolism", 
              "Replication and Repair", "Translation", "Membrane Transport", "Nucleotide Metabolism", 
              "Metabolism of Cofactors and Vitamins", "Energy Metabolism", "Amino Acid Metabolism", 
              "Carbohydrate Metabolism")
compare_genes_func$X.Datasets <- factor(compare_genes_func$X.Datasets, levels = corr_order, 
                                   labels = corr_order)

#melt your dataset
compare_genes_func_melt <- melt(compare_genes_func)
compare_genes_func_melt <-ddply(compare_genes_func_melt, .(variable), transform, rescale = scale(value))

#create palette for plot
filcols <- wes_palette("Darjeeling2", 19, type = c("continuous"))
filcols <- rainbow(4704)
#create plot with barplot, x and y labels, remove background and change theme to black lines
plot1 <- ggplot() + 
  geom_bar(aes(y = value, x = variable, fill = X.Datasets), colour = "black", data = compare_genes_func_melt, stat = "identity") +
  labs(x = "Sample Site", y = "Relative abundance (%)") +
  scale_fill_manual(values = filcols) +
  theme(axis.line.x = element_line(size=.5, colour = "black"),
        axis.line.y = element_line(size=.5, colour = "black"),
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10),
        legend.key=element_rect(fill="white", colour="white"),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())

#plot in jpeg
jpeg(filename = "barplot_genes_FUNCT_final.jpeg", width = 9000, height = 5000, res = 800)
plot1
dev.off()
