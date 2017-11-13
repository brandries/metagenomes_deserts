library(vegan)
library(ggplot2)

setwd("~/Masters/Data 2017/blast genes/Tsv files/")

otu_func <- read.delim("compare_all_genes-L3.txt", sep = "\t", header = T, row.names = 1)


otu_transf <- decostand(otu_func, "hellinger")

otu_transf <- t(otu_transf)

#bray curtis in vegan
vegdist(otu_transf, "bray") -> d

#metaMDS
fit <- metaMDS(d, "bray", k = 2, trymax = 999)

#take a look at the results
fit

#extract the scores
data.scores <- as.data.frame(scores(fit))
row3 <- c('Aus_Control_Deep','Aus_Control_Surface','Aus_FC_Deep','Aus_FC_Surface','Nam_Control_Deep_Dune','Nam_Control_Deep_GP','Nam_Control_Surface_Dune','Nam_Control_Surface_GP','Nam_FC_Deep_Dune','Nam_FC_Deep_GP','Nam_FC_Surface_Dune','Nam_FC_Surface_GP')
as.factor(row3) -> row3
row4 <- c(1,1,1,1,2,2,2,2,3,3,3,3)
row5 <- c(1,2,1,2,1,1,2,2,1,1,2,2)
row6 <- c(1,1,2,2,1,1,1,1,2,2,2,2)
as.factor(row4) -> row4
as.factor(row5) -> row5
as.factor(row6) -> row6


data.scores <- cbind(data.scores, row3,row4, row5, row6)
colnames(data.scores)[3] <- 'Site'
colnames(data.scores)[4] <- 'zone'
colnames(data.scores)[5] <- 'depth'
colnames(data.scores)[6] <- 'fc'
data.scores

#add hulls to your plot
grp.dco <- data.scores[data.scores$Site == "Dune Control", ][chull(data.scores[data.scores$Site == "Dune Control", c("NMDS1", "NMDS2")]),]
grp.dce <- data.scores[data.scores$Site == "Dune Centre", ][chull(data.scores[data.scores$Site == "Dune Centre", c("NMDS1", "NMDS2")]),]
grp.gpco <- data.scores[data.scores$Site == "Gravel Plain Control", ][chull(data.scores[data.scores$Site == "Gravel Plain Control", c("NMDS1", "NMDS2")]),]
grp.gpce <- data.scores[data.scores$Site == "Gravel Plain Centre", ][chull(data.scores[data.scores$Site == "Gravel Plain Centre", c("NMDS1", "NMDS2")]),]

hull.data <- rbind(grp.dco, grp.dce, grp.gpco, grp.gpce)
hull.data

attach(data.scores)
#plot using ggplots
plot_1 <- ggplot() +
#  geom_polygon(data = data.scores, aes(x = NMDS1, y = NMDS2, fill = depth, group = depth), alpha = 0.30) +
  scale_fill_brewer(type = div, palette = 'Set1')+
  geom_point(data = data.scores, aes(x = NMDS1, y = NMDS2, colour = Site, shape = fc), size = 6)   #add your points
  scale_color_brewer(type = div, palette = 'Set1')+
  scale_shape_manual(values = c(15, 16, 17, 18)) +
  coord_fixed(ratio = 1.05) + 
  theme_bw()  +
  guides(fill = FALSE) +
  guides(colour = FALSE) +
  guides(shape = FALSE) +
#  guides(fill = guide_legend(title = NULL))+  #remove elements of legend
#  guides(colour = guide_legend(title = NULL))+
#  guides(shape = guide_legend(title = NULL))+
  annotate("text", x = 0.23, y = -0.41, label = "stress = 0.09", size = 8) +
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


#save your plot
jpeg(filename = "bacteria2_nmds.jpeg", width = 6000, height = 5000, res = 800)
plot_nmds
dev.off()

#use multiplot
multiplot(plot_1, plot_2, plot_3, cols = 2)


######################
#THIS IS NEW, RDA
#######################

soil_chem <- read.table("./soil_chem_transformed_mgs.txt", dec = ".", row.names = 1)

rda.biol <- rda(otu_transf, soil_chem)

smry <- summary(rda.biol)
df1  <- data.frame(smry$sites[,1:2])       # PC1 and PC2

row3 <- c('Aus_Control_Deep','Aus_Control_Surface','Aus_FC_Deep','Aus_FC_Surface','Nam_Control_Deep_Dune','Nam_Control_Deep_GP','Nam_Control_Surface_Dune','Nam_Control_Surface_GP','Nam_FC_Deep_Dune','Nam_FC_Deep_GP','Nam_FC_Surface_Dune','Nam_FC_Surface_GP')
as.factor(row3) -> row3
row4 <- c(1,1,1,1,2,2,2,2,3,3,3,3)
row5 <- c(1,2,1,2,1,1,2,2,1,1,2,2)
row6 <- c(1,1,2,2,1,1,1,1,2,2,2,2)
as.factor(row4) -> row4
as.factor(row5) -> row5
as.factor(row6) -> row6
df1 <- cbind(df1, row3,row4, row5, row6)
colnames(df1)[3] <- 'Site'
colnames(df1)[4] <- 'zone'
colnames(df1)[5] <- 'depth'
colnames(df1)[6] <- 'fc'

envfit(rda.biol, soil_chem)#output gives you levels of significance
df2  <- data.frame(smry$biplot[,1:2]) 

plot_2 <- ggplot(df1, aes(x=RDA1, y=RDA2)) + 
  geom_point(aes(label=rownames(df1), colour = fc),size=6) +
  geom_segment(data=df2, aes(x=0, xend=RDA1, y=0, yend=RDA2), size = 0.5,
               color="navy", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df2, 
            aes(x=RDA1,y=RDA2,label=rownames(df2),
                hjust=0.2*(1-sign(RDA1)),vjust=0.5*(1-sign(RDA2))), 
            color="navy", size=5)
