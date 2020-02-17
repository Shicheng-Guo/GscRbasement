BiocManager::install("phyloseq")

setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/pancrease/medip")
data<-read.table("R2.txt",head=T,row.names=1,sep="\t",check.names=F)
sampleList<-unique(unlist(lapply(strsplit(colnames(data),"_"),function(x) x[1])))
library(ape)
library(phyloseq)
for(i in sampleList){
  input<-data[grep(i,colnames(data)),grep(i,colnames(data))]
  fit <- hclust(as.dist(1-input), method="ave")
  pdf(paste(i,".newick.ini.pdf",sep=""))
  par(mar=c(10,10,10,10))
  plot(fit)
  plot(dend1, nodePar = list(pch = c(1,NA), cex = 0.8, lab.cex = 0.8),type = "t", center = TRUE)
  plot(dend1, edgePar = list(col = 1:2, lty = 2:3),dLeaf = 1, edge.root = TRUE)
  plot(dend1, nodePar = list(pch = 2:1, cex = .4*2:1, col = 2:3),horiz = TRUE)
  dev.off()
  phylo <- as.phylo(fit,use.labels=TRUE)
  write.tree(phy=phylo,file=paste(i,".newick",sep=""))
}




