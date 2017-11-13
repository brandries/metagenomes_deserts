##Regression analysis for aridity##

##read in phyla and metadata file
metadata <-read.table("metadata_05March.txt",header=T, row.names=1, sep="\t", dec = ",")
bac_rel <-read.table("bac_regression.txt",header=T, row.names=1, sep="\t")
bac_rel_proteo <-read.table("bac_regression_proteo.txt",header=T, row.names=1, sep="\t")

library(vegan)
#sqrt transform to prevent biases introduced from non-normal data
decostand(bac_rel, method = "hellinger") -> bac_sqrt
decostand(bac_rel_proteo, method = "hellinger") -> bac_sqrt_proteo
write.table(bac_sqrt, file = "bac_sqrt.txt", sep = "\t", row.names = T, col.names = T)

#planc - SIG
planc<-lm(bac_sqrt$Planctomycetes~metadata$Aridity)
planc
summary(planc)
bac_sqrt$Planctomycetes~metadata$Aridity -> planc_plot

attach(metadata)
plot(planc_plot, xlab = "Aridity", ylab = "Plantomycetes", xlim = c(0.2, 1), ylim = c(0.05, 0.155), pch = 19, 
     col = c("firebrick", "aquamarine3", "dodgerblue3", "goldenrod3")[Region], 
     family="sans", font.lab=2)
detach(metadata)
abline(planc, col = "red")
svg(filename = "planctomycetes_regression.svg", height = , width = )
dev.off()

#Firmicutes - NOT SIG
firm<-lm(bac_sqrt$Firmicutes~metadata$Aridity)
firm
summary(firm)

#Cren - SIG
cren<-lm(bac_sqrt$Crenarchaeota~metadata$Aridity)
cren
summary(cren)
bac_sqrt$Crenarchaeota~metadata$Aridity -> cren_plot

attach(metadata)
plot(cren_plot, xlab = "Aridity", ylab = "Crenarchaeota", xlim = c(0.2, 1), ylim = c(0, 0.4), pch = 19, 
     col = c("firebrick", "aquamarine3", "dodgerblue3", "goldenrod3")[Region], 
     family="sans", font.lab=2)
detach(metadata)
abline(cren, col = "red")
svg(filename = "cren_regression.svg", height = , width = )
dev.off()

#Bac - SIG
bact<-lm(bac_sqrt$Bacteroidetes~metadata$Aridity)
bact
summary(bact)
bac_sqrt$Bacteroidetes~metadata$Aridity -> bac_plot

attach(metadata)
plot(bac_plot, xlab = "Aridity", ylab = "Bacteriodetes", xlim = c(0.2, 1), ylim = c(0.05, 0.45), pch = 19, 
     col = c("firebrick", "aquamarine3", "dodgerblue3", "goldenrod3")[Region], 
     family="sans", font.lab=2)
detach(metadata)
abline(bact, col = "red")
svg(filename = "bacteriodetes_regression.svg", height = , width = )
dev.off()

#Gemma - sig
gem<-lm(bac_sqrt$Gemmatimonadetes~metadata$Aridity)
gem
summary(gem)
bac_sqrt$Gemmatimonadetes~metadata$Aridity -> gem_plot

attach(metadata)
plot(gem_plot, xlab = "Aridity", ylab = "Gemmatimonadetes", xlim = c(0.2, 1), ylim = c(0.1, 0.4), pch = 19, 
     col = c("firebrick", "aquamarine3", "dodgerblue3", "goldenrod3")[Region], 
     family="sans", font.lab=2)
detach(metadata)
abline(gem, col = "red")
svg(filename = "gemm_regression.svg", height = , width = )
dev.off()

#Verr - sig
ver<-lm(bac_sqrt$Verrucomicrobia~metadata$Aridity)
ver
summary(ver)
bac_sqrt$Verrucomicrobia~metadata$Aridity -> ver_plot

attach(metadata)
plot(ver_plot, xlab = "Aridity", ylab = "Verrucomicrobia", xlim = c(0.2, 1), ylim = c(0.1, 0.45), pch = 19, 
     col = c("firebrick", "aquamarine3", "dodgerblue3", "goldenrod3")[Region], 
     family="sans", font.lab=2)
detach(metadata)
abline(ver, col = "red")
svg(filename = "verr_regression.svg", height = , width = )
dev.off()

#chloroflexi - sig
chlor<-lm(bac_sqrt$Chloroflexi~metadata$Aridity)
chlor
summary(chlor)
bac_sqrt$Chloroflexi~metadata$Aridity -> chlor_plot

attach(metadata)
plot(chlor_plot, xlab = "Aridity", ylab = "Chloroflexi", xlim = c(0.2, 1), ylim = c(0.1, 0.35), pch = 19, 
     col = c("firebrick", "aquamarine3", "dodgerblue3", "goldenrod3")[Region], 
     family="sans", font.lab=2)
detach(metadata)
abline(chlor, col = "red")
svg(filename = "chlor_regression_abline.svg", height = , width = )
dev.off()

#acidobacteria - sig
acid<-lm(bac_sqrt$Acidobacteria~metadata$Aridity)
acid
summary(acid)
bac_sqrt$Acidobacteria~metadata$Aridity -> acid_plot

attach(metadata)
plot(acid_plot, xlab = "Aridity", ylab = "Acidobacteria", xlim = c(0.2, 1), ylim = c(0.2, 0.6), pch = 19, 
     col = c("firebrick", "aquamarine3", "dodgerblue3", "goldenrod3")[Region], 
     family="sans", font.lab=2)
detach(metadata)
abline(acid, col = "red")
svg(filename = "acid_regression.svg", height = , width = )
dev.off()

#proteobacteria - not sig
proteo<-lm(bac_sqrt_proteo$Proteo~metadata$Aridity)
proteo
summary(proteo)


#actino - sig
actino<-lm(bac_sqrt$Actinobacteria~metadata$Aridity)
actino
summary(actino)
bac_sqrt$Actinobacteria~metadata$Aridity -> actino_plot

attach(metadata)
plot(actino_plot, xlab = "Aridity", ylab = "Actinobacteria", xlim = c(0.2, 1), ylim = c(0.3, 0.7), pch = 19, 
     col = c("firebrick", "aquamarine3", "dodgerblue3", "goldenrod3")[Region], 
     family="sans", font.lab=2)
detach(metadata)
abline(actino, col = "red")
svg(filename = "actino_regression.svg", height = , width = )
dev.off()

##cyano - not sig
cyan<-lm(bac_sqrt$Cyanobacteria~metadata$Aridity)
cyan
summary(cyan)

##alpha - not sig
alpha<-lm(bac_sqrt$Alphaproteobacteria~metadata$Aridity)
alpha
summary(alpha)

##beta - sig
beta<-lm(bac_sqrt$Betaproteobacteria~metadata$Aridity)
beta
summary(beta)
bac_sqrt$Betaproteobacteria~metadata$Aridity -> beta_plot

attach(metadata)
plot(beta_plot, xlab = "Aridity", ylab = "Betaproteobacteria", xlim = c(0.2, 1), ylim = c(0.1, 0.4), pch = 19, 
     col = c("firebrick", "aquamarine3", "dodgerblue3", "goldenrod3")[Region], 
     family="sans", font.lab=2)
detach(metadata)
abline(beta, col = "red")
svg(filename = "beta_regression.svg", height = , width = )
dev.off()

##delta - not sig
delta<-lm(bac_sqrt$Deltaproteobacteria~metadata$Aridity)
delta
summary(delta)
