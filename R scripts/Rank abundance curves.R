#rank abundance curves
libs <- c("ggplot2", "BiodiversityR", "plyr")
lapply(libs, require, character.only = TRUE)

#Import data
rank_data <- read.delim("~/Data 2017/blast genes/Tsv files/compare_all_tax_norm-species.txt", 
                           sep = "\t", row.names = 1)
rank_env <- read.delim("~/Data 2017/blast genes/Tsv files/compare_all_tax_norm-env.txt", 
                        sep = "\t", row.names = 1)

#transpose the frame
abundc <- t(as.matrix(rank_data))

#draw up the complete rank abundance curve
plotr <- rankabundance(abundc, y = rank_env)

#draw complete rank abundance curve
jpeg(filename = "rank_abundance_full.jpeg", width = 7000, height = 5000, res = 800)
rankabunplot(plotr,scale = "logabun", specnames=c(0:0), type = "l")
lines(plotr, lwd = 3, col = "black")

dev.off()

#draw the individual rank abundance curve
jpeg(filename = "rank_abundance_comb.jpeg", width = 7000, height = 5000, res = 800)
rankabuncomp(abundc, y = rank_env, "Sample", scale = "logabun", type = "l", legend = T)


dev.off()

svg(filename = "rank_abundance_all.svg", width = 10, height = 7)
