qqplot<-function(pvalues,output="qqplot.pdf"){
library("Haplin")
pdf(output)
  pQQ(pvalues, nlabs =length(pvalues), conf = 0.95) 
dev.off()
}
