###################
#     Q-VALUE     #  
###################

# The first step is create a txt file with p-values of our study without header. Is the same file for qq-plots.

# Load the library for extract q-values in our datasets
library("qvalue")

# Load the file with path of the file
path <- setwd("your_path_for_your_file_with_p-values.txt")

# Load the library to read input file as a table
library(data.table)

# We load our file with p-values and make a variable called "p" with that info. qobj is the comand to run our analysis with a lambda = 0
p<-fread("pvalue.txt")
qobj <-qvalue(p, lambda=0)

# In the last step, we save our results in a file in the same path as the load file.
write.csv(qobj$qvalues,"qvalue.txt")


# End of Script
