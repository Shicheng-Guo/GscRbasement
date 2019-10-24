args = commandArgs(trailingOnly=TRUE)
input=args[1]
Pcol=args[2]
output=args[3]
P=read.table(input)
pvalues=as.numeric(P[,Pcol])
qqplot<-function(pvalues,output=output){
library("Haplin")
pdf(output)
  pQQ(pvalues, nlabs =length(pvalues), conf = 0.95) 
dev.off()
}
