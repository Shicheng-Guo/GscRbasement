#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
vcfile=as.character(args[1])
output<-paste(vcfile,"pvalue.txt",sep=".")
data<-read.table(vcfile,head=T,sep="\t",skip=16,comment.char="\\",check.names=F)
phen<-read.table("~/hpc/project/pmrp/phen/IBDCH_Phetyp10_Obesity_SampleIDs.Michigen.txt",head=T)
newphen<-phen[match(colnames(data),phen[,2]),]
case<-which(newphen$PheTyp10_Obesity_C1=="2")
con<-which(newphen$PheTyp10_Obesity_C1=="1")
data<-read.table("chr22.update.vcf",head=F,sep="|",skip=17)
rlt<-c()
for(i in unique(data$INFO)){
  ca<-data[which(data$INFO %in% i),case]
  co<-data[which(data$INFO %in% i),con]
  xx<-data.matrix(restru(ca))
  yy<-data.matrix(restru(co))
  class(xx)="numeric"
  class(yy)="numeric"
  sxx<-colSums(xx)
  syy<-colSums(yy)
  A<-sum(ssum(sxx)==0)
  B<-sum(ssum(sxx)>0)
  C<-sum(ssum(syy)==0)
  D<-sum(ssum(syy)>0)
  z<-matrix(c(A,B,C,D),2,2)+0.5 
  test<-chisq.test(z)
  P=test$p.value
  temp<-data.frame(i,A,B,C,D,P,nrow(ca))
  names(temp)<-c("block","PC","NC","PN","NN","P","NSNP")
  rlt<-rbind(rlt,temp)
}
write.table(rlt,file=output,sep="\t",quote=F,col.names=T,row.names = F)
ssum<-function(vector){
rlt<-c()
  for(i in seq(1,length(vector),by=2)){
    tmp<-vector[i]+vector[i+1]
    rlt<-c(rlt,tmp)
  }
rlt
}
restru<-function(data){
  rlt<-c()
  for(i in 10:ncol(data)){
    yy<-lapply(strsplit(as.character(data[,i]),"[|]"),function(x) t(x))
    zz<-matrix(unlist(yy),ncol=2,byrow=T)
    rlt<-cbind(rlt,zz)
  }
  rlt
}
