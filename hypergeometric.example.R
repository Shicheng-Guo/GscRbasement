input<-read.table("ENSG.ENST.ENSP.Symbol.hg19.bed.txt",head=F)[,4]
P<-c()
for(j in 1:1500){
x<-unique(as.character(input))
x<-x[sample(1:length(x),4427)]
p<-sum(x %in% temp2[,1])
P<-c(P,1-phyper(p,672,24344,4427, lower.tail=T))
}
library("Haplin")
pQQ(P, nlabs =length(P), conf = 0.95)
