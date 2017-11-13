#Blast genes taxonomy analysis pipeline
#load libraries
libs <- c("wesanderson")

lapply(libs, require, character.only = TRUE)


#Import data
compare_all_tax_norm_phylum_gene <- read.delim("~/Data 2017/blast genes/Tsv files/compare_all_tax_norm-phylum_filt0.01.txt", sep = "\t", row.names = 1)
compare_all_tax_norm_phylum_read <- read.delim("~/Data 2017/blast reads/Tsv files/compare_reads_norm-phylum_filt0.01.txt", sep = "\t", row.names = 1)
compare_all_tax_norm_phylum_scaf <- read.delim("~/Data 2017/blast scaffolds/Tsv files/comp_all_scaffolds_norm-phylum-filt0.01.txt", sep = "\t", row.names = 1)

#Plot phylum gene barplot - using wes_palate

jpeg(filename = "phylum_barplot_gene.jpeg", width = 11000, height = 6000, res = 800)
barplot(t(t(compare_all_tax_norm_phylum_gene)), col =  wes_palette("Moonrise3", 12, type = c("continuous")), legend.text = TRUE, xlim = c(0,21))
dev.off()

#Plot phylum read barplot - using wes_palate

jpeg(filename = "phylum_barplot_read.jpeg", width = 11000, height = 6000, res = 800)
barplot(t(t(compare_all_tax_norm_phylum_read)), col =  wes_palette("Moonrise3", 12, type = c("continuous")), legend.text = TRUE, xlim = c(0,21))
dev.off()

#Plot phylum scaffold barplot - using wes_palate

jpeg(filename = "phylum_barplot_scaf.jpeg", width = 11000, height = 6000, res = 800)
barplot(t(t(compare_all_tax_norm_phylum_scaf)), col =  wes_palette("Moonrise3", 12, type = c("continuous")), legend.text = TRUE, xlim = c(0,21))
dev.off()

# transpose the matrix
compare_all_tax_norm_phylum_gene = t(compare_all_tax_norm_phylum_gene)

