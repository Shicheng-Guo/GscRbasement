#########################################
#     TEST OF PROPORTIONS (FISHER)      #
#########################################

# This script was buildt for R enviroment
# Remove any previous loaded objects
rm(list = ls())

# Load library to make it easier loading very large files
library(data.table)

# Define the plot settings default
opar <- par()

# Define path for files 
path <- setwd("file_with_data_for_fisher's_test")
file <- "file_with_data_for_fisher's_test"
data <- fread(file, header="auto", sep="auto", stringsAsFactors=FALSE, quote="")
head(data)
dim(data)

l <- dim(data)
l

# Create an empty variables
p.fisher <- NULL
OR.fisher <- NULL
OR.fisher.ci.left <- NULL
OR.fisher.ci.right <- NULL

for (i in 1:l) {
  
  a <- unlist(data[i,"dead_individuals_with_mutations"])
  # a
  b <- unlist(data[i,"total_dead_individuals"])
  # b
  c <- unlist(data[i,"alive_individuals_with_mutations"])
  # c
  d <- unlist(data[i,"total_alive_individuals"])
  # d
  
  # Test statistic for all the variants
  result.fisher <- fisher.test(matrix(c(a, b-a, c, d-c), ncol=2))
  result.fisher
  
  # Extract test statistic results
  p.fisher[[i]] <- result.fisher$p.value
  OR.fisher[[i]] <- result.fisher$estimate[[1]]
  OR.fisher.ci.left[[i]] <- result.fisher$conf.int[[1]]
  OR.fisher.ci.right[[i]] <- result.fisher$conf.int[[2]]
  
  # End of Loop
}

# The vector "p.fisher" have all p-values. We must consider the number of total individuals for condition 1 and for condition 2. In this case, the number is 6 and 4.
p <- format(round(p.fisher, 6), nsmall = 6)
OR <- format(round(OR.fisher, 4), nsmall = 6)
OR.ci95.left <- format(round(OR.fisher.ci.left, 6), nsmall = 6)
OR.ci95.right <- format(round(OR.fisher.ci.right, 6), nsmall = 6)

# p
# OR
# OR.ci95.left
# OR.ci95.right

# Inspect results
which(p<0.05)

# Add the p-value column to data dataframe
data2 <- cbind(data,p,OR,OR.ci95.left,OR.ci95.right)
head(data2)

# Order the results by p-value from fisher exact test of proportions
data2.ord <- data2[order(data2[,"p"]),]
head(data2.ord)
tail(data2.ord)
hist(as.numeric(data2.ord$p))

# Write results into a file to edit with EXCEL
results.to.save <- data2.ord
dim(results.to.save)

# Save the results in the same path with the name "fisher_test_of_proportions_based_on_your_variable" and complete the name file with your analysis (e.g. "fisher_test_of_proportions_based_on_mortality")
write.table(results.to.save, file="fisher_test_of_proportions_based_on_your_variable", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
