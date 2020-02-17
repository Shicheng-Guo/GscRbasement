#############################
#   GET THE VCF EXPANDED    #
#############################

# This script was buildt for R enviroment
# A vcf expanded with the number of individuals with each variant (SNV and indels) is needed. Previously, we make a recoded vcf where individuals with 0/1 or 1/0 or 1/1 are codified as "1". Individuals with 0/0 are codified as "0" and individuals without data for that variant (./.) are codified as "NA"

# Remove any previous loaded objects
rm(list = ls())

# Load library to make it easier loading very large files
library(data.table)

# Define the plot settings default 
opar <- par()

# Define path for files 
path <- setwd("your_fie_path")
file <- "input_vcf_recode.vcf"
data <- fread(file, header="auto", sep="auto", stringsAsFactors=FALSE, quote="")
head(data)
dim(data)

# Define the affected and unaffected individuals in your study.
# Condition1: Variable1 is "1" in samples X1, X2, X3, X7 and X9.
# Condition2: Variable1 is "0" in samples X4, X5, X6, X8 and X10. 

# Number of sequenced individuals
# Update these values according to the covariate under study. Data for "condition1" are (for that example):
ind.cond1 <- 5
ind.cond2 <- 5

# Sum the number of individuals in condition1 having a "1" genotype (mutant genotype corresponding to 0/1 or 1/0 or 1/1)
sum.measured.geno.cond1 <- rowSums(data[,c("X1","X2","X3","X7","X9")], na.rm=TRUE)
sum.measured.geno.cond1
length(sum.measured.geno.cond1)

# Get the number of "NA" in columns corresponding to died individuals (condition "1")
sum.na.geno.cond1 <- rowSums(is.na(data[,("X1","X2","X3","X7","X9")]))
sum.na.geno.cond1
length(sum.na.geno.cond1)

# The number of usable condition 1 individuals at each row is (ind.cond1)-(number of NA's in individuals with this condition)
sum.usable.geno.cond1 <- ind.cond1 - sum.na.geno.cond1
sum.usable.geno.cond1
length(sum.usable.geno.cond1)

# Sum the number of condition2 individuals having a "0" genotype (non-mutant genotype or 0/0)
sum.measured.geno.cond2 <- rowSums(data[,c("X4"","X5"","X6"","X8","X10")], na.rm=TRUE)
sum.measured.geno.cond2
length(sum.measured.geno.cond2)

# Get the number of "NA" in columns corresponding to alive individuals (condition "2")
sum.na.geno.cond2 <- rowSums(is.na(data[,c("X4"","X5"","X6"","X8","X10")]))
sum.na.geno.cond2
length(sum.na.geno.cond2)

# The number if usable alive (condition "2") individuals at each row is (ind.cond2)-(number of NA's in individuals with this condition)
sum.usable.geno.cond2 <- ind.cond2 - sum.na.geno.cond2
sum.usable.geno.cond2
length(sum.usable.geno.cond2)

# Add these computed variables to the original data matrix into a number of new colums
# "measured" stands for sequenced genotypes
# "usable" stands for usable genotypes after removal of "NA" genotypes (absent genotypes)
data2 <- cbind(data,sum.measured.geno.cond1,sum.usable.geno.cond1,sum.measured.geno.cond2,sum.usable.geno.cond2)
head(data2)
dim(data2)
# Our variant data is contained in "data2"


# The second file that we need for this analysis is a gene list with the chr, start position, end position and gene name (HUGO name). The order of the columns must be as is indicated previously. This list can be dowload in UCSC brownser. The order of the columns in that file must be as is described in
# Load detailed genes list
path <- setwd("your_gene_list_variants_path")
file <- "your_genes_with_variants_details.txt"
genes <- fread(file, header="auto", sep="auto", stringsAsFactors=FALSE, quote="")
head(genes)
dim(genes)

# Get the number of genes
l <- dim(genes)[1]
l

# Get the number of variants (be aware of a possible header)
m <- dim(data2)[1]
m

# Create an empty vector
genename <- NULL

for (i in 1:l) {
  genes.chr <- genes[i,1]
  genes.start <- genes[i,2]
  genes.stop <- genes[i,3]
  genes.name <- genes[i,4]
  
  # genes.chr
  # genes.start
  # genes.stop
  # genes.name
  
  message(date(),"; Gene iteration:",i,"; looking for variants at ",genes.chr," within gene ",genes.name)
  log <- paste(date(),"; Gene iteration:",i,"; working at ",genes.chr," with gene ",genes.name)
  write(log, file = "log_Gene_Enrichment_Analysis_Rcode_1", append = TRUE)
  
  # Iterave over the variants in the VCF and count genotypes in each condition
  # Then annotate the gene (if the variant is within a known gene)
  
  for (j in 1:m) {
    
    variant.chr <- data2[j,1]
    variant.chr
    variant.pos <- data2[j,2]
    variant.pos
    
    # message("Variant iteration: ",j)
    
    # Then, compare genes.chr and variant.chr, and genes.start and variant.pos
    # If they match, then annotate the gene into this variant
    # Else proceed to the next variant until the end of the variant list
    # Iterate with the next gene in the genes list
    
    if (variant.chr == genes.chr) {
      
      if ( (variant.pos >= genes.start) && (variant.pos <= genes.stop) ) {
        message("Found a variant within gene ",genes.name)
        # annotate gene name into the variant
        genename[j] <- genes.name
      }
    }
    
    # End of FOR loop for variants: do not delete it! 
  }
  
  # End of FOR loop for genes: do not delete it! 
}

# Check length of column-vector "genename"
length(genename)

# Check if there is any NULL value
length(unlist(genename))

# Make a copy of genename
gn <- genename
length(gn)
length(unlist(gn))
gn[sapply(gn, is.null)] <- "nogene"
length(unlist(gn))

# Add the column-vector "genename" (it contains the genes corresponding to the variants) to the data2 matrix
# Remember to unlist the list of genes to convert it into a column-vector
data3 <- cbind(data2,unlist(gn))

# Set colname for "genename" column
names(data3)[24]<-"genename"

# Save data into an expanded VCF
write.table(data3, file="somatics_recode_expanded.vcf", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
