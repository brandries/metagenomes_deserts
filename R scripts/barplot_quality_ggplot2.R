#Barplots for quality filtering etc
libs <- c("ggplot2", "ggthemes", "extrafont", "plyr", "scales", "wesanderson")
lapply(libs, require, character.only = TRUE)

#load data for barplots
barplot_data <- read.delim("C:/Users/Andries van der Walt/Dropbox/Masters 2017/Thesis/Data and analysis FC/barplots_seq_quality_f.txt", 
                           sep = "\t")

#set levels and arrange factors in the correct order for the graph as well as rename the legend
barplot_data$Type <- factor(barplot_data$Type, levels = c("a - Input sequences", "b - Good sequences","c - BLAST hits" ), 
                            labels = c("Bad Sequences", "Good Sequencing", "BLAST hits found"))

corr_location <- c("ACD", "ACS",
                   "NCDD",
                   "NCSD",
                   "NCDG",
                   "NCSG",
                   "AFD",
                   "AFS",
                   "NFDD",
                   "NFSD",
                   "NFDG" ,
                   "NFSG"
)
compare_genes_tax$Location<- factor(compare_genes_tax$Location, levels = corr_location, labels = corr_location)

#create palette for plot
filcols <- wes_palette("Darjeeling2", 3)

#create plot with barplot, x and y labels, remove background and change theme to black lines
plot1 <- ggplot() + 
  geom_bar(aes(y = Sequences, x = Location, fill = Type), colour = "black", data = barplot_data, stat = "identity") +
  labs(x = "Sample Site", y = "Sequence number") +
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
jpeg(filename = "barplot_quality_final.jpeg", width = 7000, height = 5000, res = 800)
plot1
dev.off()

#plot in svg
svg(filename = "barplot_quality_fc.svg", width = 8, height = 4, pointsize = 12)
plot1
dev.off()
