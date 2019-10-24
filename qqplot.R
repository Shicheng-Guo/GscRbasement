args = commandArgs(trailingOnly=TRUE)
library("Haplin")
input=args[1]
Pcol=as.numeric(args[2])
data=read.table(input,head=T,sep="")
pvalues=as.numeric(data[,Pcol])
qqplot<-function(pvalues,output=paste(input,"qqplot.pdf",sep="")){
pdf(output)
  pQQ(pvalues, nlabs =length(pvalues), conf = 0.95)
dev.off()
}
# wget https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/qqplot.R -O qqplot.R
# Rscript manhattan.plot.R Schrodi_IL23_IL17.Kbac.assoc 4 Schrodi_IL23_IL17.Kbac.assoc
