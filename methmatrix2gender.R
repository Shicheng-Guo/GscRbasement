#!/usr/bin/env Rscript

###  Prediction gender with methylation matrix
# rowname is coordination
# colnames is sample id
# First version MHB regions don't have any overlap with FHM therefore can not make gender inference.
# 2017-05-23
######################################################################################################
##########################load functions ##############################################################
######################################################################################################
#install.packages("Deducer")
#install.packages("stringr")
#install.packages("pROC")
library("stringr")
library("pROC")
#library("Deducer")

cor2bed<-function(cor){
  cor<-as.character(cor)
  a<-unlist(lapply(strsplit(cor,split=c(":")),function(x) strsplit(x,"-")))
  bed<-matrix(a,ncol=3,byrow=T)
  bed<-data.frame(bed,cor)
  return(data.frame(bed))
}
bed2cg<-function(bed1){
  ref<-read.table("~/oasis/db/hg19/GPL13534.sort.bed",head=F,sep="\t")
  cor2bed<-function(cor){
    cor<-as.character(cor)
    a<-unlist(lapply(strsplit(cor,split=c(":")),function(x) strsplit(x,"-")))
    bed<-matrix(a,ncol=3,byrow=T)
    bed<-data.frame(bed,cor)
    return(data.frame(bed))
  }
  rbedintersect<-function(bed1,ref){
    Rbedtools<-function(functionstring="intersectBed",bed1,bed2,opt.string=""){
      #create temp files
      a.file=tempfile()
      b.file=tempfile()
      out   =tempfile()
      options(scipen =99) # not to use scientific notation when writing out
      #write bed formatted dataframes to tempfile
      write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
      write.table(bed2,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
      # create the command string and call the command using system()
      command=paste(functionstring,"-a",a.file,"-b",b.file,opt.string,">",out,sep=" ")
      cat(command,"\n")
      try(system(command))
      res=read.table(out,header=F)
      unlink(a.file);unlink(b.file);unlink(out)
      res=subset(res,V5!=".")
      return(res)
    }
    merge<-Rbedtools(functionstring="intersectBed",bed1,ref,opt.string="-wao")
    return(merge)
  }
  merge<-rbedintersect(bed1,ref)
  return(merge)
}

bed2cor<-function(bed,extend=0){
  cor<-apply(bed,1,function(x) paste(x[1],":",as.numeric(x[2])-extend,"-",as.numeric(x[3])+extend,sep=""))
  cor<-gsub(" ","",cor)
  return(cor)
}

cg2bed<-function(cg,extend=0){
  bed2cor<-function(bed){
    cor<-apply(bed,1,function(x) paste(x[1],":",as.numeric(x[2])-extend,"-",as.numeric(x[3])+extend,sep=""))
    cor<-gsub(" ","",cor)
    return(cor)
  }
  ref<-read.table("~/oasis/db/hg19/GPL13534.sort.bed",head=F,sep="\t")
  bed<-ref[match(cg,ref[,4]),1:3]
  bed[,2]=bed[,2]-extend
  bed[,3]=bed[,3]+extend
  cor<-bed2cor(bed)
  rlt<-data.frame(bed,cor,cg)
  return(rlt)
}

write.bed<-function(bed,file,extend=0){
  bed[,2]<-as.numeric(as.character(bed[,2]))-extend
  bed[,3]<-as.numeric(as.character(bed[,3]))+extend
  if(ncol(bed)==3){
    bed[,4]<-paste(bed[,1],":",bed[,2],"-",bed[,3],sep="")  
  }
  if(ncol(bed)>=4){
    write.table(bed,file=file,sep="\t",col.names=F,row.names=F,quote=F)
  }
}

see<-function(x){
  x[1:3,1:3]
}

readmeth450<-function(){
  rlt<-list()
  library("stringr")
  file<-list.files(pattern="jhu*")
  data<-c()
  for(i in file){
    tmp<-read.table(i,head=T,skip=1,row.names=1,sep="\t",check.names = FALSE,as.is=T)
    data<-cbind(data,tmp[,1])
    print(i)
  }
  colnames(data)<-unlist(lapply(unlist(lapply(file,function(x) unlist(strsplit(x,"[.]"))[6])),function(x) substr(x,1,15)))
  rownames(data)<-rownames(tmp)
  cancertype<-unique(unlist(lapply(file,function(x) unlist(strsplit(x,"_|.Human"))[2])))
  sampletype<-unlist(lapply(unlist(lapply(file,function(x) unlist(strsplit(x,"[.]"))[6])),function(x) substr(x,14,15)))
  save(data,file=paste(cancertype,"meth.RData",sep="."))
  rlt$data<-data
  rlt$cancertype<-cancertype
  rlt$sampletype<-sampletype
  rlt$cpg<-rownames(data)
  return(rlt)
}

rbedintersect<-function(bed1,ref){
  Rbedtools<-function(functionstring="intersectBed",bed1,bed2,opt.string=""){
    #create temp files
    a.file=tempfile()
    b.file=tempfile()
    out   =tempfile()
    options(scipen =99) # not to use scientific notation when writing out
    #write bed formatted dataframes to tempfile
    write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
    write.table(bed2,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
    # create the command string and call the command using system()
    command=paste(functionstring,"-a",a.file,"-b",b.file,opt.string,">",out,sep=" ")
    cat(command,"\n")
    try(system(command))
    res=read.table(out,header=F)
    unlink(a.file);unlink(b.file);unlink(out)
    res=subset(res,V5!=".")
    return(res)
  }
  merge<-Rbedtools(functionstring="intersectBed",bed1,ref,opt.string="-wao")
  return(merge)
}
methmatrix2gender<-function(data){
  probability<-apply(data,2,function(x) sum(na.omit(x)>0.3)/length(na.omit(x)))
  gender<-c(NA,length(probability))
  gender[probability<0.3]<-"male"
  gender[probability>0.3]<-"female"
  prediction<-data.frame(id=substr(as.character(names(probability)),1,12),gender,probability)
  return(prediction)
}
#########################################################################################
########################## prediction ####################################
#########################################################################################
args = commandArgs(trailingOnly=TRUE)
methmatrix<-args[1]
data<-read.table(methmatrix,head=T,row.names=1,sep="\t",as.is=T)
fhm<-read.table("~/oasis/db/FHM.bed",head=F,sep="\t",as.is=T)
testdata<-data[match(unique(bed2cor(rbedintersect(cor2bed(rownames(data)),fhm)[,1:3])),rownames(data)),]
predict<-methmatrix2gender(testdata)
print(predict)

