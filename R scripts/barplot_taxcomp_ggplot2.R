#Barplots for taxonomic composition
libs <- c("ggplot2", "ggthemes", "extrafont", "plyr", "scales", "wesanderson")
lapply(libs, require, character.only = TRUE)

#Import data
#Your data has to be a table with first column being the one of axis, and the second being the other axis.
#the third column contains the values used to draw the bar chart. 
#alternatively, you can input a OTU table style table, and change the format using the reshape::melt function.

compare_genes_tax <- read.delim("./compare_all_tax_norm-phylum_filt0.01_r.txt", sep = "\t", row.names = 1)

#order the figure - pur in the phyla in your own dataset
corr_order <- c("Other", "Ascomycota", "Acidobacteria", "Actinobacteria ", "Bacteroidetes", "Chloroflexi", "Cyanobacteria",
                "Firmicutes", "Gemmatimonadetes", "Planctomycetes", "Proteobacteria", "Thaumarchaeota")
corr_nam <- c("Other", "Ascomycota", "Acidobacteria", "Actinobacteria", "Bacteroidetes", "Chloroflexi", "Cyanobacteria",
                "Firmicutes", "Gemmatimonadetes", "Planctomycetes", "Proteobacteria", "Thaumarchaeota")
corr_location <- c("ACD", "ACS",
                   "AFD",
                   "AFS",
                   "NCDD",
                   "NCSD",
                   "NCDG",
                   "NCSG",
                   "NFDD",
                   "NFSD",
                   "NFDG" ,
                   "NFSG"
)
compare_genes_tax$Phylum <- factor(compare_genes_tax$Phylum, levels = corr_order, 
                            labels = corr_nam)
compare_genes_tax$Location<- factor(compare_genes_tax$Location, levels = corr_location)
#create palette for plot
filcols <- wes_palette("Darjeeling2", 12, type = c("continuous"))
filcols <- rainbow(12)

#create plot with barplot, x and y labels, remove background and change theme to black lines
plot1 <- ggplot() + 
  geom_bar(aes(y = Abundance, x = Location, fill = Phylum), colour = "black", data = compare_genes_tax, stat = "identity") +
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
        panel.background = element_blank())+
  facet_grid(.~Source , scales = "free")

#plot in jpeg
jpeg(filename = "barplot_genes_tax_facet.jpeg", width = 7000, height = 5000, res = 800)
plot1
dev.off()

#plot in svg
svg(filename = "barplot_genes_facet_fc.svg", width = 8, height = 4, pointsize = 12)
plot1
dev.off()
