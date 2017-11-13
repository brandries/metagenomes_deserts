#   THIS IS A SHORT SCRIPT TO QUERY STRINGS IN ONE DATA FRAME AGAINST A LARGER DATABASE AND RETREIVE 
#   CORRESPONDING VALUES FROM A DIFFERENT COLUMN IN THE DATAFRAME
#   THIS EXAMPLE EXTRACTS SIGNIFICANT COOCURRING OTUS FROM A LARGER SET OF OTUS

#IMPORT DATA INTO R ENVIORNMENT
#TABLE IS THE DATABASE
#NAME IS THE QUERY NAMES

setwd("C:/Users/Andries van der Walt/Dropbox/Masters 2017/Thesis/CHapter 2 - FC interactions/AMplicon analysis/otu_tables")


otu_database <-read.table("combined.namib.16S.fasta_tax_assignments.txt",header=F, sep="\t")
otu_query <-read.table("1515F.fastq.trimmed.filt.trim.fasta.biom.txt",header=F, sep="\t")

#SET COUNTING VARIABLES TO 0

x = 0
y = 0

#CREATE EMPTY MATRIX WITH TWO COLOUMNS AND AS MANY ROWS AS THERE ARE QUERIES
matrix(nrow = length(otu_database$V1), ncol = 31) -> otu_query_taxonomy

#EMBEDDED FOR LOOPS WITH AN IF STATEMENT SEQUENTIALLY TESTING IF THE QUERY STRING IS IN THE DATABASE

for(i in otu_database$V1){
  y = y+1
  x = 0
  for(a in otu_query$V1){
    x = x+1
    if(as.character(i) == as.character(otu_query$V1[x])){
      otu_query_taxonomy[y,1] <- i
      otu_query_taxonomy[y,2] <- as.character(otu_query$V2[x])
      
    } 
    
  }
}

#EXPORT TABLE

write.table(otu_query_taxonomy, file = "otu_query_taxonomy.txt", sep = "\t")

