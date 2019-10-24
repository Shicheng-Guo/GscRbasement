args = commandArgs(trailingOnly=TRUE)
input=args[1]
Pcol=args[2]
P=read.table(input)
pvalues=as.numeric(P[,Pcol])
qqplot<-function(pvalues,output=paste(input,"qqplot.pdf",sep=""){
library("Haplin")
pdf(output)
  pQQ(pvalues, nlabs =length(pvalues), conf = 0.95) 
dev.off()
}
# wget https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/qqplot.R -O qqplot.R
# Rscript manhattan.plot.R Schrodi_IL23_IL17.Kbac.assoc 4 Schrodi_IL23_IL17.Kbac.assoc

