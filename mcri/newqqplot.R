logcminput<-function(memo){
  memo="RA500"
  library("CMplot")
  file<-list.files(pattern="*logistic.adjusted")
  output<-read.table(file[1],head=T)
  chr2num<-function(x){
    x<-output$V1
    x<-gsub("chr","",x)
    x[x=="X"]<-23
    x[x=="Y"]<-24
    return(x)
  }
  POS<-unlist(lapply(strsplit(as.character(output$SNP),":"),function(x) x[2]))
  cminput<-na.omit(data.frame(SNP=output$SNP,Chromosome=output$CHR,Position=POS,trait1=output$UNADJ,stringsAsFactors=F))
  CMplot(cminput,plot.type="b",memo="logistic",LOG10=TRUE,threshold=NULL,file="jpg",dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6)
  write.csv(cminput,file=paste(memo,"logistic.cminput.txt",sep=""),quote=F,row.names=F)
}

assocminput<-function(memo){
  output<-read.table("RA500.counts.assoc",head=T)
  chr2num<-function(x){
    x<-output$V1
    x<-gsub("chr","",x)
    x[x=="X"]<-23
    x[x=="Y"]<-24
    return(x)
  }
  cminput<-na.omit(data.frame(SNP=output$SNP,Chromosome=output$CHR,Position=output$BP,trait1=output$P,stringsAsFactors=F))
  CMplot(cminput,plot.type="b",memo="assoc",LOG10=TRUE,threshold=NULL,file="jpg",dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6)
  write.csv(cminput,file=paste("assoc.cminput.txt",sep=""),quote=F,row.names=F)
}


