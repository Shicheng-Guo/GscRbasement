#####################
#     QQ-PLOTS      #
#####################

# Extract in bash the column with p-values to perform qq-plots. In this case, the comand in bash is:
awk '{print $9}' fisher_test_of_proportions_based_on_your_variable > pvalue_for_your_variable.txt

# Then, we load in R the p-values list. 
path <- setwd("path_of_your_pvalue_for_your_variable.txt")
file <- "pvalue_for_your_variable.txt"

pvals <- read.table("pvalue_mortality.txt", header=T)

# It is needed to know the name of the column for p-values (in this case, the name is "p").
observed <- sort(pvals$p)
lobs <- -(log10(observed))

expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))

# A PDF is created with the qq-plot in the same path.
pdf("qqplot.pdf", width=6, height=6)
plot(c(0,7), c(0,7), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,7), ylim=c(0,7), las=1, xaxs="i", yaxs="i", bty="l")
points(lexp, lobs, pch=23, cex=.4, bg="black") 
dev.off()
