## Prediction gender with 450K dataset
# Download json and transfer to csv (https://json-csv.com/)
# Only can be run in Linux since it requre bedtools

#install.packages("Deducer")
install.packages("stringr")
install.packages("pROC")
library("stringr")
library("pROC")
#library("Deducer")

bed2cg<-function(bed1){
  ref<-read.table("~/work/db/hg19/GPL13534.sort.bed",head=F,sep="\t")
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

cg2bed<-function(cg){
  bed2cor<-function(bed){
    cor<-apply(bed,1,function(x) paste(x[1],":",x[2],"-",x[3],sep=""))
    cor<-gsub(" ","",cor)
    return(cor)
  }
  ref<-read.table("~/work/db/hg19/GPL13534.sort.bed",head=F,sep="\t")
  bed<-ref[match(cg,ref[,4]),1:3]
  cor<-bed2cor(bed)
  rlt<-data.frame(bed,cor)
  return(rlt)
}

setwd("/media/NAS3_volume2/shg047/HM450/TCGA/lihc")
saminfo<-read.csv("../clinical.project-TCGA-LIHC.2017-05-23T04-50-49.306129.csv")
saminfo<-data.frame(id=substr(as.character(saminfo$exposures__submitter_id),1,12),gender=as.character(saminfo$demographic__gender))

library("stringr")
file<-list.files(pattern="jhu*")
data<-c()
for(i in file){
  tmp<-read.table(i,head=T,skip=1,row.names=1,sep="\t",check.names = FALSE,as.is=T)
  data<-cbind(data,tmp[,1])
  print(i)
}

#load("PancancerMethMatrix_March2016.RData")
#load("PancancerMethMatrix_March2016.Test.RData")
colnames(data)<-unlist(lapply(unlist(lapply(file,function(x) unlist(strsplit(x,"[.]"))[6])),function(x) substr(x,1,15)))
rownames(data)<-rownames(tmp)
cancertype<-unique(unlist(lapply(file,function(x) unlist(strsplit(x,"_|.Human"))[2])))
sampletype<-unlist(lapply(unlist(lapply(file,function(x) unlist(strsplit(x,"[.]"))[6])),function(x) substr(x,14,15)))
save(data,file=paste(cancertype,"meth.RData",sep="."))
# idv<-as.array(str_extract(file,"TCGA-[0-9|a-z|A-Z]*-[0-9|a-z|A-Z]*-[0-9]*"))
# length(idv)
# colnames(data)<-idv
bed<-cg2bed(rownames(data))
newdata1<-data[which(bed[,1]=="chrX"),]
newdata2<-data[which(bed[,1]=="chrY"),]
dim(newdata1)
dim(newdata2)
gender<-saminfo[match(substr(colnames(data),1,12),saminfo[,1]),2]
delta<-data.frame(t(apply(newdata1,1,function(x) tapply(x,gender,function(x) mean(x,na.rm=T)))))
marker<-rownames(subset(delta,female<0.6 & female>0.4 & male<0.1))
test<-data[match(marker,rownames(data)),]
input<-data.frame(score=apply(test,2,function(x) sum(x>0.3,na.rm=T)/(length(na.omit(x)))),gender)
fit<-glm(gender~score,input,family=binomial(link = "logit"))
prob=predict(fit,type=c("response"))
input$prob=prob
g <- roc(gender ~ prob, data = input)
jpeg("ROC1.jpg")
plot(g)  

# give up since 
#modelfit <- glm(formula=gender ~ score, family=binomial(), data=input, na.action=na.omit)
#jpeg("ROC2.jpg")
#rocplot(modelfit)
#dev.off()
