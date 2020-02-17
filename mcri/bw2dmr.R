#!/usr/bin/env Rscript
setwd("/mnt/bigdata/Genetic/Projects/shg047/db/hg19/wgbs/slideWindow/")
TtestPValue<-function(data,x1,x2,pair=FALSE){
  data<-data.matrix(data)
  output<-matrix(NA,dim(data)[1],6)   # set output matrix ()
  for(i in 1:dim(data)[1]){
    out<-data.frame()
    if(pair==TRUE){
      Valid<-nrow(na.omit(data.frame(data[i,x1],data[i,x2])))
    }else{
      Valid<-100
    }
    if( sum(!is.na(data[i,x1]))>=3 & sum(!is.na(data[i,x2]))>=3 & Valid>3){ 
      tmp1<-try(t.test((data[i,x1]),(data[i,x2]),paired=F, na.action=na.omit))
      output[i,1]<-format(tmp1$p.value, scientific=TRUE)
      output[i,2]<-round(mean((data[i,x1]))-mean((data[i,x2])),3)
      output[i,3]<-round(mean((data[i,x1])),3)
      output[i,4]<-round(mean((data[i,x2])),3)
      output[i,5]<-round(sd(data[i,x1]),3)
      output[i,6]<-round(sd(data[i,x2]),3)
      print(i)
    }
  }
  rownames(output)<-rownames(data)
  output
}

args = commandArgs(trailingOnly=TRUE)
sam<-read.table("/gpfs/home/guosa/hpc/epimarker/cell.index.txt",head=F,sep="\t")
blood=which(sam[,4]=="Blood")
solid=which(sam[,4]=="Solid")
rlt<-TtestPValue(args[1],blood,solid)
output=paste(args[1],".pvalue.rlt",sep="")
write.table(rlt,file=output,sep="\t",quote=F,col.names=T,row.names=F)


