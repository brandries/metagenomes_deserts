#this is a script for drawing a scatterplot and then regressing over it
libs <- c("ggplot2", "plyr", "wesanderson", "RColorBrewer")
lapply(libs, require, character.only = TRUE)

setwd("C:/Users/Andries van der Walt/Dropbox/Masters 2017/Thesis/Data and analysis FC/")
quality <- read.delim("Metagenomes_quality.txt", row.names = 1, header = T)

#create palette for plot
filcols <- wes_palette("Darjeeling2", 2, type = c("continuous"))
##THESE WERE ALL TESTS

  #create your model
  modl <- lm(quality$Alignment_Rate_...~quality$Length_.bp.1000)
  modl
  summary(modl)

  xmin <- min(quality$Predicted_Genes)
  xmax <- max(quality$Predicted_Genes)
  predicted <- data.frame(Predicted_Genes=seq(xmin, xmax, length.out = 12))

  predicted$Predicted_Genes_assigned_taxonomy <- predict(modl, predicted)
  predicted



  #another test
  sp <- ggplot(quality, aes(x=Predicted_Genes, y=Predicted_Genes_assigned_taxonomy, colour=Predicted_Genes, shape=Location, fill=Depth))+
    geom_point(size=4)+                          #this sets the point size and shape unless you set shape manually
    scale_shape_manual(values = c(21,22,23)) +        #this sets the shape values
    scale_color_continuous(low = "blue", high = "red")+
    scale_fill_manual(values = c("orange", "purple"))+
    geom_smooth(method = lm)
  
  sp  
  sp + geom_line(data=predicted, size = 1) + geom_point(size = 2)



  #TEST with basics
  sp <- ggplot(quality, aes(x=Predicted_Genes, y=Predicted_Genes_assigned_taxonomy)) +
    geom_point(colour="violet")
  
  sp + geom_line(data=predicted, size = 1)

####################################
#THIS IS WHERE THE REAL SCRIPT STARTS
####################################
  
  
#THIS IS THE REAL PLOT genes vs assembly span

attach(quality)
modl <- lm(quality$Predicted_Genes~quality$Alignment_Rate_...)
modl
summary(modl)

sp <- ggplot(quality, aes(x=Length_.bp.1000, y=Predicted_Genes, fill= "orange2"))+
  geom_point(shape=21, size=4)+                          #this sets the point size and shape unless you set shape manually
  scale_color_manual(values = "black")+
  scale_fill_manual(values = "orange2")+
  guides(fill=FALSE)+
  scale_x_continuous(limits = c(0,60000000))+
  geom_smooth(method = lm, se = F, colour = "royalblue4")+
  annotate("text", x=60000000, y=-0.05, label = "R^2 == -0.095", parse=T) +
  annotate("text", x=60000000, y=-50000, label = "p == 0.83", parse=T)+
  labs(x="% Short reads aligned", y = "Number of genes predicted")+
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

sp


#safe the plot
svg(filename = "span_vs_genes.svg", width = 6, height = 4, pointsize = 12)

sp

dev.off()

#THIS IS THE REAL PLOT for alignment vs contigs
attach(quality)
modl <- lm(quality$Alignment_Rate_...~quality$X._contigs_.1000)
modl
summary(modl)

sp <- ggplot(quality, aes(x=X._contigs_.1000, y=Alignment_Rate_..., fill= "orange2"))+
  geom_point(shape=21, size=4)+                          #this sets the point size and shape unless you set shape manually
  scale_color_manual(values = "black")+
  scale_fill_manual(values = "orange2")+
  guides(fill=FALSE)+
  scale_x_continuous(limits = c(0,30000))+
  geom_smooth(method = lm, se = F, colour = "royalblue4")+
  annotate("text", x=30000, y=-0.05, label = "R^2 == 0.603", parse=T) +
  annotate("text", x=30000, y=-5, label = "p == 0.002", parse=T)+
  labs(x="Number of contigs generated", y = "% Short reads aligned")+
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

sp


#safe the plot
svg(filename = "alignment_v_contigs.svg", width = 6, height = 4, pointsize = 12)

sp

dev.off()

#THIS IS THE REAL PLOT for alignment vs genes preditced
attach(quality)
modl <- lm(quality$Predicted_Genes~quality$Alignment_Rate_...)
modl
summary(modl)

sp <- ggplot(quality, aes(x=Alignment_Rate_..., y=Predicted_Genes, fill= "orange2"))+
  geom_point(shape=21, size=4)+                          #this sets the point size and shape unless you set shape manually
  scale_color_manual(values = "black")+
  scale_fill_manual(values = "orange2")+
  guides(fill=FALSE)+
  scale_x_continuous(limits = c(0,60))+
  geom_smooth(method = lm, se = F, colour = "royalblue4")+
  annotate("text", x=60, y=-0.05, label = "R^2 == -0.095", parse=T) +
  annotate("text", x=60, y=-50000, label = "p == 0.830", parse=T)+
  labs(y="Number of predicted genes", x = "% Short reads aligned")+
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

sp


#safe the plot
svg(filename = "alignment_v_genes.svg", width = 6, height = 4, pointsize = 12)

sp

dev.off()


#THIS IS THE REAL PLOT contigs vs genes preditced
attach(quality)
modl <- lm(quality$Predicted_Genes~quality$X._contigs_.1000)
modl
summary(modl)

sp <- ggplot(quality, aes(x=X._contigs_.1000, y=Predicted_Genes, fill= "orange2"))+
  geom_point(shape=21, size=4)+                          #this sets the point size and shape unless you set shape manually
  scale_color_manual(values = "black")+
  scale_fill_manual(values = "orange2")+
  guides(fill=FALSE)+
  scale_x_continuous(limits = c(0,30000))+
  geom_smooth(method = lm, se = F, colour = "royalblue4")+
  annotate("text", x=30000, y=-0.05, label = "R^2 == 0.205", parse=T) +
  annotate("text", x=30000, y=-50000, label = "p == 0.078", parse=T)+
  labs(y="Number of predicted genes", x = "Number of contigs generated")+
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

sp


#safe the plot
svg(filename = "contigs_v_genes.svg", width = 6, height = 4, pointsize = 12)

sp

dev.off()
