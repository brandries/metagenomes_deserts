#Script to get otu table into phyloseq

libs <- c("ggplot2", "seqtime", "SpiecEasi", "plyr",  "reshape", "Matrix", "phyloseq", "ape", "igraph", "BIOMASS")
lapply(libs, require, character.only = TRUE)


setwd("C:/Users/Andries van der Walt/Google Drive/Masters/Data 2017/blast genes/Tsv files/")

#read the otu table, mapping file and reference taxonomy sets

tax <- read.delim("./compare_all_genes_norm-L3_order_nohuman.txt", sep = "\t", header = T, row.names = 1)
mapping_file <- read.delim("./metadata_andreis.txt", row.names = 1)
taxonomy_mapping <- read.delim("./mapping_function.txt", row.names = 1)

#transpose the otu table and make the reference taxa a matrix
taxonomy_mapping <- as.matrix(taxonomy_mapping)

#create otu table and tax for phyloseq objects and make the combined object this is for the whole community
OTU <- otu_table(otu_namib, taxa_are_rows = T)
TAX <- tax_table(taxonomy_mapping)
SAMPL <- sample_data(mapping_file)
physeq = phyloseq(OTU, TAX, SAMPL)
physeq

otu_namib <-tax[,5:12]
#network
spiec.out=spiec.easi(physeq, method="mb", lambda.min.ratio = 1e-2, nlambda = 20,icov.select.params=list(rep.num=20))
spiec.graph=adj2igraph(spiec.out$refit, vertex.attr=list(name=taxa_names(physeq)))
plot_network(spiec.graph, physeq, type='taxa', color="L1", label=NULL)


#split the networks for the FC and Control samples
otu_control <- tax[,1:6]
otu_fc <-tax[,7:12]

#create otu table and tax for phyloseq objects and make the combined object this is for the control community
OTU <- otu_table(otu_control, taxa_are_rows = T)
TAX <- tax_table(taxonomy_mapping)
SAMPL <- sample_data(mapping_file)
physeq = phyloseq(OTU, TAX, SAMPL)
physeq

#network
spiec.out=spiec.easi(physeq, method="mb", lambda.min.ratio = 1e-2, nlambda = 20,icov.select.params=list(rep.num=20))
spiec.graph=adj2igraph(spiec.out$refit, vertex.attr=list(name=taxa_names(physeq)))
plot_network(spiec.graph, physeq, type='taxa', color="L2", label=NULL)

#create otu table and tax for phyloseq objects and make the combined object this is for the control community
OTU <- otu_table(otu_fc, taxa_are_rows = T)
TAX <- tax_table(taxonomy_mapping)
SAMPL <- sample_data(mapping_file)
physeq = phyloseq(OTU, TAX, SAMPL)
physeq

#network
spiec.out=spiec.easi(physeq, method="mb", lambda.min.ratio = 1e-2, nlambda = 20,icov.select.params=list(rep.num=20))
spiec.graph=adj2igraph(spiec.out$refit, vertex.attr=list(name=taxa_names(physeq)))
plot_network(spiec.graph, physeq, type='taxa', color="L1", label=NULL)


#Export to cytoscape

write.graph(spiec.graph,file="spieceasi.ncol_namib.txt",format="ncol") 
write.table(TAX,file="taxonomy_network.txt",sep="\t", quote=FALSE)

?#get positive and negative interaction
betaMat=as.matrix(symBeta(getOptBeta(spiec.out)))


otu.ids=colnames(spiec.out$data)
edges=E(spiec.graph)
edge.colors=c()
for(e.index in 1:length(edges)){
  adj.nodes=ends(spiec.graph,edges[e.index])
  xindex=which(otu.ids==adj.nodes[1])
  yindex=which(otu.ids==adj.nodes[2])
  beta=betaMat[xindex,yindex]
  if(beta>0){
    edge.colors=append(edge.colors,"forestgreen")
  }else if(beta<0){
    edge.colors=append(edge.colors,"red")
  }
}
E(spiec.graph)$color=edge.colors

#now plot using red and green for positive and negative
#this is not very effective and should not be used
spiec.graph.b=spiec.graph
nodenames=V(spiec.graph.b)$name
V(spiec.graph.b)$name=getTaxonomy(nodenames, taxa.f, useRownames=TRUE)
E(spiec.graph.b)$arrow.size=5
V(spiec.graph.b)$color="white"
V(spiec.graph.b)$frame.color="black"
tkplot(spiec.graph.b)

#cluster the network
clusters=cluster_fast_greedy(spiec.graph)
clusterOneIndices=which(clusters$membership==1)
clusterOneOtus=clusters$names[clusterOneIndices]
clusterTwoIndices=which(clusters$membership==2)
clusterTwoOtus=clusters$names[clusterTwoIndices]
clusterThreeIndices=which(clusters$membership==3)
clusterThreeOtus=clusters$names[clusterThreeIndices]
clusterFourIndices=which(clusters$membership==4)
clusterFourOtus=clusters$names[clusterFourIndices]
clusterFiveIndices=which(clusters$membership==5)
clusterFiveOtus=clusters$names[clusterFiveIndices]

sort(table(getTaxonomy(clusterOneOtus,taxa.f,useRownames = TRUE)),decreasing = TRUE)

#Export to cytoscape
write.graph(spiec.graph,file="spieceasi.ncol.txt",format="ncol") 
write.table(TAX,file="taxonomy_network.txt",sep="\t", quote=FALSE)
