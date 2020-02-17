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
  manhattan(manhattaninput,col = c("blue4", "orange3"),ylim = c(0,9),genomewideline=-log10(genomewideline),lwd=1.5, suggestiveline=F)
  dev.off()
}



miRdata<-read.table("/home/guosa/hpc/db/mirTarget/Predicted_Target_Locations.default_predictions.hg19.bed")
RAGene<-read.table("/home/guosa/hpc/rheumatology/RA/KEGG_90_RA_GeneList.txt")


miR<-miRdata[miRdata[,4] %in% RAGene[,1],]
write.table(miR,file="miRNA_KEGG_90Gene.txt",sep="\t",quote=F,col.names = F,row.names=F)

setwd("/gpfs/home/guosa/hpc/rheumatology/RA/miRNA")
d1<-read.table("KEGG_miRNA_Target.txt")
d2<-read.table("RA_miRNA_SNP.txt")
d3<-read.table("KEGG_RA_miRNA.txt")

d1[d1[,5] %in% d2[,1],]
d2[d2[,1] %in% d3[,1],]


write.table(d1[d1[,5] %in% d3[,1],],file="miRNA_RA_KEGG.hg10.txt",sep="\t",quote=F,col.names = F,row.names=F)
