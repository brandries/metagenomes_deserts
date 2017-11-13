libs <- c("ggplot2", "vegan", "plyr",  "reshape" )
lapply(libs, require, character.only = TRUE)

setwd("C:/Users/Andries van der Walt/Google Drive/Masters/Data 2017/MAPLE metabolic potential")

maple_file <- read.delim("./MAPLE_ALL_pathways_R.txt", dec = ".", row.names = 1)
mapping_file <- read.delim("./mapping.txt")

attach(mapping_file)

adonis(maple_file~Location)


heatmap(as.matrix(t(maple_file[,1:140])))
heatmapr(maple_file)


library(superheat)
maple_file_t <- as.data.frame(t(maple_file))
maple_file <- as.numeric(maple_file[3:14,])



svg(filename = "heatmap_maple.svg", width = 24, height = 12, pointsize = 10)
superheat(maple_file,  pretty.order.rows = T, clustering.method = "kmeans",smooth.heat = TRUE)
dev.off()

#plot in pdf  For some reason svg does not work
pdf(file = "heatmap_maple.pdf")

#plot in jpeg
png(filename = "heatmap_maple.png", width = 7000, height = 6000, res = 800)

