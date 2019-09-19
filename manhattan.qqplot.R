args = commandArgs(trailingOnly=TRUE)
mylimma<-read.table(args,head=T,sep="",check.names=F)
seed=sample(seq(1,100000,by=1),1)
manhattan.plot<-function(mylimma){
library(qqman)
CHR=mylimma$CHR
if(length(grep("X",CHR))>0){
  CHR<-sapply(CHR,function(x) gsub(pattern = "X",replacement = "23",x))
  CHR<-sapply(CHR,function(x) gsub(pattern = "Y",replacement = "24",x))
}
CHR<-as.numeric(CHR)
manhattaninput=data.frame(SNP=mylimma$SNP,CHR=CHR,BP=mylimma$BP,P=mylimma$P)
max<-max(2-log(manhattaninput$P,10))
genomewideline=0.05/nrow(manhattaninput)
pdf(paste("manhattan.",seed,".pdf",sep=""))
manhattan(manhattaninput,col = c("blue4", "orange3"),ylim = c(0,10),genomewideline=-log10(genomewideline),lwd=2, suggestiveline=F)
dev.off()
}
manhattan.plot(mylimma)

qqplot<-function(pvalues,output="qqplot.pdf"){
library("Haplin")
pdf(paste("qqplot.",seed,".pdf",sep=""))
  pQQ(na.omit(pvalues), nlabs =length(pvalues), conf = 0.99)
dev.off()
}
qqplot(mylimma$P)

# usage
# plink --bfile RA2020-B8 --pca --threads 31
# plink --bfile RA2020-B8 --logistic --covar plink.eigenvec --covar-number 1-5 --adjust
# grep "ADD\|NMISS" plink.assoc.logistic > plink.assoc.logistic.add
# wget https://raw.githubusercontent.com/Shicheng-Guo/ASA/master/manhattan.plot.R -O manhattan.plot.R
# Rscript manhattan.plot.R plink.assoc.logistic.add
