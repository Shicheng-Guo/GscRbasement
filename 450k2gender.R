#########################################################################################
##########################  Shicheng Guo #############################################
#########################################################################################
## Prediction gender with 450K dataset
# Download json and transfer to csv (https://json-csv.com/)
# Only can be run in Linux since it requre bedtools
#########################################################################################
##########################Training with LIHC Dataset ####################################
#########################################################################################
#install.packages("Deducer")
#install.packages("stringr")
#install.packages("pROC")
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

cg2bed<-function(cg,extend=100){
  bed2cor<-function(bed){
    cor<-apply(bed,1,function(x) paste(x[1],":",as.numeric(x[2])-extend,"-",as.numeric(x[3])+extend,sep=""))
    cor<-gsub(" ","",cor)
    return(cor)
  }
  ref<-read.table("~/work/db/hg19/GPL13534.sort.bed",head=F,sep="\t")
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
  #load("PancancerMethMatrix_March2016.RData")
  #load("PancancerMethMatrix_March2016.Test.RData")
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


setwd("/media/NAS3_volume2/shg047/HM450/TCGA/lihc")
saminfo<-read.csv("../clinical.project-TCGA-LIHC.2017-05-23T04-50-49.306129.csv")
saminfo<-data.frame(id=substr(as.character(saminfo$exposures__submitter_id),1,12),gender=as.character(saminfo$demographic__gender))



load("LIHC.meth.RData")
# idv<-as.array(str_extract(file,"TCGA-[0-9|a-z|A-Z]*-[0-9|a-z|A-Z]*-[0-9]*"))
# length(idv)
# colnames(data)<-idv
bed<-cg2bed(rownames(data))
newdata1<-data[which(bed[,1]=="chrX"),]
dim(newdata1)
gender<-saminfo[match(substr(colnames(data),1,12),saminfo[,1]),2]
delta<-data.frame(t(apply(newdata1,1,function(x) tapply(x,gender,function(x) mean(x,na.rm=T)))))
marker<-rownames(subset(delta,female<0.6 & female>0.4 & male<0.1))
test<-data[match(marker,rownames(data)),]
input<-data.frame(femaleScore=apply(test,2,function(x) sum(x>0.3,na.rm=T)/(length(na.omit(x)))),gender)

bed<-cg2bed(marker)
 
F1<-newdata1[,grep("female",gender)]
M1<-newdata1[,grep("male",gender)]
png("density.png")
plot(density(na.omit(as.numeric(M1))),col="blue",lwd=2)
lines(density(na.omit(as.numeric(F1))),col="red",lwd=2)
legend("topright",legend=c("female","male"),lwd=2,col=c("red","blue"))
dev.off()

## evaluation
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
#########################################################################################
##########################Test with ESCA Dataset ####################################
#########################################################################################
setwd("/media/NAS3_volume2/shg047/HM450/TCGA/chol")
saminfo<-read.csv("../clinical.project-TCGA-CHOL.2017-05-23T06-36-48.852956.csv")
saminfo<-data.frame(id=substr(as.character(saminfo$exposures__submitter_id),1,12),gender=as.character(saminfo$demographic__gender))
chol<-readmeth450()
data<-chol$data
gender<-data.frame(id=substr(as.character(colnames(data)),1,12),gender=unlist(saminfo[match(substr(colnames(data),1,12),saminfo[,1]),2]))

testdata2gender<-function(data){
  probability<-apply(data[match(bed$cg,rownames(data)),],2,function(x) sum(na.omit(x)>0.3)/length(na.omit(x)))
  gender<-c(0,length(probability))
  gender[probability<0.3]<-"male"
  gender[probability>0.3]<-"female"
  prediction<-data.frame(id=substr(as.character(names(probability)),1,12),gender,probability)
  return(prediction)
}
predict<-testdata2gender(data)
rlt<-merge(predict,gender,by="id")
rlt[which(! as.character(rlt$gender.x)==as.character(rlt$gender.y)),]
# only one sample (TCGA-W5-AA2T) prediction error, but I do think the sample was labelled to female, however, the truth of that sample is male samples. 


