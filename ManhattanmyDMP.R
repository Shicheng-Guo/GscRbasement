ManhattanmyDMP<-function(myDMP){
  library(qqman)
  SNP=rownames(myDMP)
  CHR=myDMP$CHR
  if(length(grep("X",CHR))>0){
    CHR<-sapply(CHR,function(x) gsub(pattern = "X",replacement = "23",x))
    CHR<-sapply(CHR,function(x) gsub(pattern = "Y",replacement = "24",x))
  }
  CHR<-as.numeric(CHR)
  BP=myDMP$MAPINFO
  P=myDMP$P.Value
  manhattaninput=data.frame(SNP,CHR,BP,P)
  max<-max(2-log(manhattaninput$P,10))
  genomewideline=0.05/nrow(manhattaninput)
  genomewideline=max(subset(myDMP,adj.P.Val<0.05)$P.Value)
  pdf("manhattan.pdf")
  manhattan(manhattaninput,col = c("blue4", "orange3"),ylim = c(0,10),lwd=2, suggestiveline=F)
  dev.off()
}
