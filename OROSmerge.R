OROSmerge<-function(memo){
memo="housekeeping"
file<-list.files(pattern="*pvalue.csv")
dge<-read.csv(file[1])
os<-read.csv(file[2])
out<-merge(dge,os,by="X")

pick<-subset(out,pval<10^-4 & pval.fixed<0.05 & beta*(TE.fixed-1)>0)
pick<-pick[order(pick$pval),]

if(nrow(pick)>80){
pick<-pick[1:80,]
}

dir.create("pick")
# be careful to run it or else it will takes quite long time
for(i in 1:nrow(pick)){
x<-paste("cp *-",pick[i,1],"* ./pick",sep="")
system(x)
}
# install.packages("CMplot")
setwd("pick")
system("cp ../*.csv ./")
library("CMplot")
ensg<-read.table("~/hpc/db/hg19/ENSG.hg19.bed")
ensgrlt<-ensg[match(unlist(lapply(strsplit(as.character(out$X),"[.]"),function(x) x[1])),ensg[,5]),]
output<-data.frame(ensgrlt,out)

chr2num<-function(x){
x<-output$V1
x<-gsub("chr","",x)
x[x=="X"]<-23
x[x=="Y"]<-24
return(x)
}

write.csv(out,file=paste(memo,".tcga.pancancer.overall.rnaseq.dmg.smd",".os.hr.csv",sep=""),quote=F)
write.csv(pick,file=paste(memo,".tcga.pancancer.pick.rnaseq.dmg.smd","os.hr.csv",sep=""),quote=F)

cminput<-na.omit(data.frame(SNP=output[,ncol(output)],Chromosome=chr2num(output$V1),Position=output$V2,trait1=output$pval,stringsAsFactors=F))
CMplot(cminput,plot.type="b",memo=paste(memo,".dge.fix",sep=""),LOG10=TRUE,threshold=NULL,file="jpg",dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6)
write.csv(cminput[order(cminput[,ncol(cminput)]),],file=paste(memo,".pval.manhattan.qqplot.meta.dge.bed",sep=""),quote=F,row.names=F)

cminput<-na.omit(data.frame(SNP=output[,ncol(output)],Chromosome=chr2num(output$V1),Position=output$V2,trait1=output$pval.fix))
CMplot(cminput,plot.type="b",memo=paste(memo,".fix",sep=""),LOG10=TRUE,threshold=NULL,file="jpg",dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6)
write.csv(cminput[order(cminput[,ncol(cminput)]),],file=paste(memo,".pval.fix.manhattan.qqplot.meta.OS.HR.bed",sep=""),quote=F,row.names=F)

cminput<-na.omit(data.frame(SNP=output[,ncol(output)],Chromosome=chr2num(output$V1),Position=output$V2,trait1=output$pval.random))
CMplot(cminput,plot.type="b",memo=paste(memo,".random",sep=""),LOG10=TRUE,threshold=NULL,file="jpg",dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6)
write.csv(cminput[order(cminput[,ncol(cminput)]),],file=paste(memo,".pval.random.manhattan.qqplot.meta.OS.HR.bed",sep=""),quote=F,row.names=F)
}
