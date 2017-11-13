libs <- c("ggplot2", "ggthemes", "BiodiversityR", "vegan", "reshape")
lapply(libs, require, character.only = TRUE)


#Import data from functional
setwd("C:/Users/Andries van der Walt/Google Drive/Masters/Data 2017/blast genes/Tsv files/")
otu_func <- read.delim("./compare_all_genes-L2_nohuman.txt", sep = "\t", header = T, row.names = 1)
otu_func <- read.delim("compare_all_genes_norm-L3_order_nohuman.txt", sep = "\t", header = T, row.names = 1)
#Import soilchem
setwd("C:/Users/Andries van der Walt/Google Drive/Masters/Data 2017/Soil Chemistry")
soil_chem <- read.delim("./soil_chem_final_all_reorder.txt", dec = ".", row.names = 1)
mapping_file <- read.delim("./metadata_mgs.txt")
#import taxa
setwd("C:/Users/Andries van der Walt/Google Drive/Masters/Data 2017/blast genes/Tsv files/")
otu_tax <- read.delim("./compare_all_tax_norm-phylum_filt0.01.txt", sep = "\t", header = T, row.names = 1)



soil_dist <- vegdist(soil_chem, method = "euclidean")
pca1 = prcomp(soil_chem, scale. = TRUE)    # Set principal components

attach(mapping_file)
scores = as.data.frame(pca1$x)
summary(pca1)                 #gives cumulative variance explained

# plot of observations
shps = c(22, 21, 24)  

plot_pca <- ggplot() +
  geom_point(data = scores, aes(x = PC1, y = PC2, label=rownames(scores), fill = Loc_dep, shape = Source), size = 6) +  #add your points
  scale_fill_brewer(type = "div", palette = "Paired")+
  scale_shape_manual(values = shps) +
  coord_fixed(ratio = 1.85) + 
  theme_bw()  +
  xlab("PC1 (58.67%)") +
  ylab("PC2 (19.92%)") +
  theme(axis.text.x = element_blank(), #most of arguments to follow is to remove features of the graph
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 18), #increase the size of x label text
        axis.title.y = element_text(size = 18), #increase size of y label text
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18))

#Test significance
adonis(soil_chem~ Location)



#save your plot
svg(filename = "soil_chem_pca_mgs.svg", width = 12, height = 6, pointsize = 12)
plot_pca
dev.off()

##########################################################
#Compare to fuctional data
##########################################################
library(Hmisc)
library(corrplot)
#library(ggcorrplot)
#library(viridis)

otu_func <- decostand(otu_func, method = "hellinger")
otu_tax <- decostand(otu_tax, method = "hellinger")
soil_chem <- decostand(soil_chem, method = "hellinger")


corrtest <- as.matrix(cbind(soil_chem, t(otu_tax), t(otu_func)))

#correlation_reslut <- cor(corrtest)
#correlation_p <- cor_pmat(corrtest)

rcorr_result <- rcorr(corrtest)

svg(filename = "correlations_wtax.svg", width = 12, height = 12, pointsize = 14)
corrplot(rcorr_result$r, type = "upper", tl.col="black", tl.cex = 0.7, method = "circle",
         p.mat = rcorr_result$P, insig = "blank", diag=FALSE, sig.level = 0.05)
dev.off()

jpeg(filename = "correlations.jpeg", width = 10000, height = 5000, res = 800)
