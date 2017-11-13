#Barplots for taxonomic composition
libs <- c("ggplot2", "wesanderson", "reshape", "dplyr", "RColorBrewer", "viridis", "vegan", "plyr")
lapply(libs, require, character.only = TRUE)

#Import data
setwd("C:/Users/Andries van der Walt/Google Drive/Masters/Data 2017/blast genes/Tsv files/")
kegg_heatmap <- read.delim("./compare_all_genes_norm-L3.txt", sep = "\t", header = T, row.names = 1)

kegg_heatmap <- decostand(kegg_heatmap, "log")
#reorder dataframe, setting X.datasets as factor and order by first column 
kegg_heatmap$X.Datasets <- with(kegg_heatmap, reorder(X.Datasets, Nam_Control_Deep_Dune))

#change dimentions into three rows, with values in third column and then transform and rescale data
kegg_heatmap_m <- melt(kegg_heatmap)
kegg_heatmap_m <- ddply(kegg_heatmap_m, .(variable), transform, rescale = scale(value))

plotheatmap <- ggplot(kegg_heatmap_m, aes(variable, X.Datasets)) + geom_tile(aes(fill = rescale), colour = "white") + 
  scale_fill_gradient(low = "blue", high = "red")
#different palette
 scale_fill_viridis(name="# Events", label=comma)

#OR
  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
  scale_fill_gradient(colours = myPalette(11))

jpeg(filename = "heatmap.jpeg", width = 7000, height = 5000, res = 800)
plotheatmap
dev.off()
