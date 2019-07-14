ManhattanPlot<-function(mylimma){
library(qqman)
res <- mylimma
SNP=res$probeID
CHR=res$CHR
if(length(grep("X",CHR))>0){
  CHR<-sapply(CHR,function(x) gsub(pattern = "X",replacement = "23",x))
  CHR<-sapply(CHR,function(x) gsub(pattern = "Y",replacement = "24",x))
}
CHR<-as.numeric(CHR)
BP=res$MAPINFO
P=res$P.Value
manhattaninput=data.frame(SNP,CHR,BP,P)
max<-max(2-log(manhattaninput$P,10))
genomewideline=0.05/nrow(manhattaninput)
pdf("manhattan.pdf")
manhattan(manhattaninput,col = c("blue4", "orange3"),ylim = c(0,10),genomewideline=-log10(genomewideline),lwd=2, suggestiveline=F)
dev.off()
}
