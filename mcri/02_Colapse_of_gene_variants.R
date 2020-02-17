#######################################
#   COLAPSING THE GENE VARIANTS 
#######################################

# This script was buildt for R enviroment
# Remove any previous loaded objects
rm(list = ls())

# Load library to make it easier loading very large files
library(data.table)

# Define the plot settings default
opar <- par()

# Define path for files 
path <- setwd("path_of_your_vcf_expanded")
file <- "your_vcf_input_recode_expanded.vcf"
data <- fread(file, header="auto", sep="auto", stringsAsFactors=FALSE, quote="")
head(data)
dim(data)

# Remove any previous results
system("rm gene_enrichment_results")

# Define initial values
sum.cond1 <- 0
sum.ind.cond1 <- 0
sum.cond2 <- 0
sum.ind.cond2 <- 0

l <- dim(data)[1]
l

# Create an empty list
# dato <- NULL
gn <- NULL

for (i in 1:l) {
  
  chr <- data[i,"#CHROM"]
  genename <- data[i,"genename"]
  
  # Extract data
  sum.cond1 <- sum.cond1 + data[i,"sum.measured.geno.cond1"]
  sum.ind.cond1 <- sum.ind.cond1 + data[i,"sum.usable.geno.cond1"]
  sum.cond2 <- sum.cond2 + data[i,"sum.measured.geno.cond2"]
  sum.ind.cond2 <- sum.ind.cond2 + data[i,"sum.usable.geno.cond2"]
  
  message("i= ",i)
  
  if ( i < l ) {
    
    if ( data[i,"genename"] != data[i+1,"genename"] && data[i,"genename"] != "nogene" ) {
      message("Next genename is distinct than actual...")
      p1 <- format(round(sum.cond1/sum.ind.cond1, 4), nsmall = 4)
      p2 <- format(round(sum.cond2/sum.ind.cond2, 4), nsmall = 4)
      dato <- c(chr,genename,sum.cond1,sum.ind.cond1,sum.cond2,sum.ind.cond2,p1,p2)
      gn[[i]] <- genename
      write.table(dato, file="gene_enrichment_results", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE, append=TRUE)
      
      # Reset the sums
      sum.cond1 <- 0
      sum.ind.cond1 <- 0
      sum.cond2 <- 0
      sum.ind.cond2 <- 0
    }
  }
  
  if ( i == l && data[i,"genename"] != "nogene" ) {
    p1 <- format(round(sum.cond1/sum.ind.cond1, 4), nsmall = 4)
    p2 <- format(round(sum.cond2/sum.ind.cond2, 4), nsmall = 4)
    dato <- c(chr,genename,sum.cond1,sum.ind.cond1,sum.cond2,sum.ind.cond2,p1,p2)
    write.table(dato, file="gene_enrichment_results", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE, append=TRUE)
    gn[[i]] <- genename
  }
  
  # End of FOR loop
}

header <- paste("chr","genename","sum.ind.mutant.cond1","sum.ind.cond1","sum.ind.mutant.cond2","sum.ind.cond2","proportion.mutant.cond1","proportion.mutant.cond2",sep="\t")
write.table(header, file="header",row.names=FALSE, col.names=FALSE, quote=FALSE)
system("cat header gene_enrichment_results > gene_enrichment_results_fnl")
