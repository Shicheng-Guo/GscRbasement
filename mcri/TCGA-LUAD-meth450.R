source("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/Pancancer_mh450/meth450Pancancer.R")
source("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/pancancer/methylation/meth450Pancancer.R")
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/GscTools.R")
##################################################################################################### 
####################### Step 1: Read methylation 450K files ######################################### 
##################################################################################################### 
setwd("/mnt/bigdata/Genetic/Projects/shg047/methylation/Pancancer")
load("methdata.pancancer.RData")
methdata[1:5,1:5]
phen4<-id2phen4(colnames(methdata))
phen3<-id2phen3(colnames(methdata))
bin<-id2bin(colnames(methdata))
pid<-id2pid(colnames(methdata))
N<-c()
iid<-unique(pid)
iid<-iid[-14]
for(i in iid){
  print(i)
  input<-methdata[,which(pid==i)]
  bin<-id2bin(colnames(input))
  input<-input[,c(which(bin==1),which(bin==11))]
  print(dim(input))
  bin<-id2bin(colnames(input))
  N<-rbind(N,table(bin))
}
rownames(N)<-iid
##################################################################################################### 
####################### Step 2: Select Cancer Type ################################################## 
#####################################################################################################
library("SIS")
library("arm")
library("randomForest")

DMR<-function(methdata,pid=""){
  pid<-id2pid(colnames(methdata))
  input<-methdata[,which(pid=="LUAD"|pid=="LUSC")]
  phen4<-id2phen4(colnames(input))
  phen3<-id2phen3(colnames(input))
  bin<-id2bin(colnames(input))
  pid<-id2pid(colnames(input))
  input<-input[,c(which(bin==1),which(bin==11))]
  phen4<-id2phen4(colnames(input))
  phen3<-id2phen3(colnames(input))
  bin<-id2bin(colnames(input))
  pid<-id2pid(colnames(input))
  phen<-data.frame(phen4=phen4,phen3=phen3,pid=pid,bin=bin)
  head(phen)
  map<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/hg19/GPL13534_450K_hg19_V3.bed")
  system("wget https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/Pancancer_mh450/Normal.PBMC.GEO.HM450K.Beta.RData")
  load("Normal.PBMC.GEO.HM450K.Beta.RData")
  system("wc -l ~/hpc/db/hg19/GPL13534_450K_hg19_PBMC_BUR.bed")
  BUR<-read.table("~/hpc/db/hg19/GPL13534_450K_hg19_PBMC_BUR.bed")
  input<-input[rownames(input) %in% BUR$V4,]
  PDMR<-read.table("~/hpc/methylation/Pancancer/TCGA-Pancancer-MH450.Meta.diff.txt",head=T,row.names=1,sep="\t")
  DMR<-subset(PDMR,beta>0.1 & pval<10^-5)
  DMG<-na.omit(cpg2symbol(rownames(DMR)))
  N<-length(unique(DMG))
  input<-input[rownames(input)%in%rownames(DMR),]
  bin<-id2bin(colnames(input))
  input<-data.frame(phen=bin,t(input))
  input<-data.frame(t(na.omit(t(input))))
  input$phen[input$phen==11]<-0
  P=apply(input[,2:ncol(input)],2,function(x) summary(bayesglm(as.factor(input[,1])~x,family=binomial))$coefficients[2,4])
  newinput<-input[,c(1,match(names(P[head(order(P),n=6000)]),colnames(input)))]
  save(newinput,file="newinput.RData")
  RF <- randomForest(as.factor(phen) ~ ., data=newinput, ntree=5000,importance=TRUE,proximity=T)
  imp<-RF$importance
  imp<-imp[order(imp[,4],decreasing = T),]
  newimp<-data.frame(imp,Symbol=cpg2symbol(row.names(imp)))
  newimp<-na.omit(newimp)
  RLT<-data.frame(newimp,map[match(rownames(newimp),map[,4]),])
  START=RLT$V2-150
  END=RLT$V3+150
  CHR=RLT$V1
  RLT<-data.frame(CHR,START,END,RLT)
  write.table(RLT,file=paste("../../RF_LUAD_LUSC.BUR.PAN.Top6000VIP.txt",sep=""),sep="\t",quote=F,col.names = NA,row.names = T)
}


for(i in 1:10){
TEMP<-newinput[,c(1,match("cg05336395",colnames(newinput)))]
RF <- randomForest(as.factor(phen) ~ ., data=TEMP, ntree=5000,importance=TRUE,proximity=T)
}

##################################################################################################### 
####################### Step 4: Tumor Suppressor GENES ############################################### 
##################################################################################################### 
setwd("~/hpc/methylation/Pancancer/RNA-seq")
library("arm")
library("randomForest")
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/GscTools.R")
load("rnaseqdata.pancancer.RData")
phen1=read.table("~/hpc/methylation/TCGA-clinical-11093.tsv",header = T,sep="\t")
phen2=read.table("~/hpc/methylation/Pancancer/RNA-seq/File_metadata2.txt",header = T,sep="\t")
barcode=read.table("~/hpc/project/TCGA/pancancer/FPKM/barcode.txt",head=T,sep="\t")
phen<-data.frame(phen1[match(phen2$cases.0.case_id,phen1$case_id),],phen2)
phen$file_name=gsub(".gz","",phen$file_name)
colnames(rnaseqdata)<-unlist(lapply(strsplit(colnames(rnaseqdata),"/"),function(x) x[2]))
phen<-phen[match(colnames(rnaseqdata),phen$file_name),]
phen$phen4<-id2phen4(phen$cases.0.samples.0.submitter_id)
phen$phen3<-id2phen3(phen$cases.0.samples.0.submitter_id)
phen$bin<-id2bin(phen$cases.0.samples.0.submitter_id)
phen$pid=phen$project_id

for( i in c("KIRC","BRCA","MESO","THCA","HNSC","PRAD","LIHC","UCEC","KIRP","LUSC","COAD","LUAD","BLCA","ESCA","PAAD","OV","CHOL","READ","SARC")){
  ppid=paste("TCGA-",i,sep="")
  print(ppid)
  rlt<-DGE(rnaseqdata,pid=ppid)
}


DGE<-function(rnaseqdata,pid="TCGA-LUSC"){
  yphen<-subset(phen,(bin==1 | bin==11) & pid==pid)
  input<-rnaseqdata[,match(yphen$file_name,colnames(rnaseqdata))]
  dim(yphen)
  dim(input)
  colnames(input)<-yphen$phen4
  Seq<-paste(yphen$pid,yphen$bin,sep="-")
  y<-abs(as.numeric(as.factor(yphen$bin))-2)
  rlt<-c()
  for(i in 1:nrow(input)){
    fit<-summary(bayesglm(y~input[i,]))$coefficients[2,]
    Mean<-tapply(input[i,],y,function(x) mean(x,na.rm=T))
    fit<-c(Mean,fit)
    rlt<-rbind(rlt,fit)
  }
  rownames(rlt)<-rownames(input)
  write.table(rlt,file=paste("~/hpc/methylation/",pid,"-RNAseq-FPKM-UQ.ENSG.DEG.txt",sep=""),sep="\t",col.names = T,row.names=T,quote=F)
  rownames(rlt)<-ENSG2Symbol(rownames(input))
  write.table(rlt,file=paste("~/hpc/methylation/",pid,"-RNAseq-FPKM-UQ.Symbol.DEG.txt",sep=""),sep="\t",col.names = T,row.names=T,quote=F)
  save(rlt,file=paste("~/hpc/methylation/",pid,"-RNAseq-FPKM-UQ.Symbol.DEG.RData",sep=""))
  return(rlt)
}

temp<-data.frame(RLT,rlt[match(RLT$Symbol,rownames(rlt)),])
output<-subset(temp,t.value<0 & Pr...t..<10^-5)
write.table(output,file=paste("RF.",LUAD,".BUR.PAN.Top4000VIP.txt",sep=""),col.names = NA,row.names = T,quote=F)
##################################################################################################### 
####################### Step 4: Tumor Suppressor GENES ############################################### 
##################################################################################################### 
setwd("/mnt/bigdata/Genetic/Projects/shg047/methylation/Pancancer")
newimp<-head(RLT,n=200)
TSG<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/TSGene2.0.txt",head=T,sep="\t")
TSGTARGET<-newimp[newimp$Symbol %in% TSG[,2],]
write.table(TSGTARGET,file=paste("../../RF_LUAD_LUSC.BUR.PAN.Top6000VIP.TSG.txt",sep=""),sep="\t",quote=F,col.names = NA,row.names = T)

newimp<-head(RLT,n=170)
GWAS<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/lung/lung_GWAS.Genelist.txt")
DMRGWAS<-na.omit(data.frame(newimp,GWAS[match(newimp$Symbol,GWAS[,1]),]))
unique(DMRGWAS$Symbol)
write.table(DMRGWAS,file=paste("../../RF_LUAD_LUSC.BUR.PAN.Top6000VIP.WGBS.txt",sep=""),sep="\t",quote=F,col.names = NA,row.names = T)

newimp<-head(RLT,n=3000)
DEG<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/Pancancer_FPKM_UQ/TCGA-Pancancer-RNAseq-FPKM-UQ.Meta.diff.Symbol.txt",sep="\t",head=T)
DMRDEG<-data.frame(newimp,DEG[match(newimp$Symbol,DEG$hgnc_symbol),])
DMRDEGTarget<-subset(DMRDEG,beta<0 & pval<10^-5)
unique(DMRDEGTarget$hgnc_symbol)
ADOS<-read.table("~/hpc/methylation/TCGA_LUAD_FPKM-UQ.DGE_OS_HR_PanDiff.All.txt",head=T,sep="\t")
SCOS<-read.table("~/hpc/methylation/TCGA_LUSC_FPKM-UQ.DGE_OS_HR_PanDiff.All.txt",head=T,sep="\t")
ADOSS<-subset(ADOS,Estimate<0 & beta<0 & Pr...t..<10^-5 & Pr...z..<0.01)
SCOSS<-subset(SCOS,Estimate<0 & beta<0 & Pr...t..<10^-5 & Pr...z..<0.01)
xx<-newimp[which(newimp$Symbol %in% ADOSS$Symbol),]
yy<-newimp[which(newimp$Symbol %in% SCOSS$Symbol),]
MethPANExpOS<-unique(rbind(xx,yy))
unique(MethPANExpOS$Symbol)
write.table(MethPANExpOS,file=paste("../../RF_LUAD_LUSC.BUR.PAN.Top6000VIP.OS.txt",sep=""),sep="\t",quote=F,col.names = NA,row.names = T)

##################################################################################################### 
####################### Step 5: Tumor Suppressor GENES ############################################### 
##################################################################################################### 
set.seed(49)
cv.error <- NULL
k <- 10
rlt1<-c()
rlt2<-c()
for(i in 1:k){
  index <- sample(1:nrow(input),round(0.9*nrow(input)))
  train.cv <- input[index,]
  test.cv <- input[-index,]
  P=apply(train.cv[,2:ncol(train.cv)],2,function(x) summary(bayesglm(as.factor(train.cv[,1])~x,family=binomial))$coefficients[2,4])
  train.cv<-train.cv[,c(1,match(names(P[head(order(P),n=200)]),colnames(train.cv)))]
  test.cv<-test.cv[,c(1,match(names(P[head(order(P),n=200)]),colnames(test.cv)))]
  meth<-list()
  meth$train.cv=train.cv
  meth$test.cv=test.cv
  print(paste(ncol(train.cv),"variables passed P-value threshold and enrolled in SIS model"))
  RF <- randomForest(as.factor(phen) ~ ., data=train.cv, importance=TRUE,proximity=T)
  imp<-RF$importance
  imp<-imp[order(imp[,4],decreasing = T),]
  head(imp)
  write.table(imp,file=paste("RandomForest.VIP.Meth.",i,".txt",sep=""),sep="\t",quote=F,row.names = T,col.names = NA)
  topvar<-match(rownames(imp)[1:30],colnames(input))
  train.cv <- input[index,c(1,topvar)]
  test.cv <- input[-index,c(1,topvar)]
  n <- colnames(train.cv)
  f <- as.formula(paste("phen ~", paste(n[!n %in% "phen"], collapse = " + ")))
  nn <- neuralnet(f,data=train.cv,hidden=c(5,3),act.fct = "logistic",linear.output = F,threshold = 0.1)
  pr.nn <- neuralnet::compute(nn,test.cv)
  trainRlt<-data.frame(phen=train.cv[,1],pred=unlist(nn$net.result[[1]][,1]))
  testRlt<-data.frame(phen=test.cv[,1],pred=unlist(pr.nn$net.result[,1]))
  rownames(trainRlt)=row.names(train.cv)
  rownames(testRlt)=row.names(test.cv)
  rlt1<-rbind(rlt1,trainRlt)  
  rlt2<-rbind(rlt2,testRlt)
  print(i)
}
data1<-na.omit(data.frame(rlt1))
data2<-na.omit(data.frame(rlt2))
model.glm1 <- bayesglm(phen~.,data=rlt1,family=binomial(),na.action=na.omit)
model.glm2 <- bayesglm(phen~.,data=rlt2,family=binomial(),na.action=na.omit)
pred1 <- predRisk(model.glm1)
pred2 <- predRisk(model.glm2)
par(mfrow=c(2,2),cex.lab=1.5,cex.axis=1.5)
plotROC(data=data1,cOutcome=1,predrisk=cbind(pred1))
plotROC(data=data2,cOutcome=1,predrisk=cbind(pred2))
##################################################################################################### 
####################### Step 6: TSNE to show subtype ############################################### 
##################################################################################################### 
tsneplot(input[,match(rownames(RLT),colnames(input)),2:ncol(input)],input[,1],plot="TCGA.BRCA.BUR.TSNE.TopVariable.pdf",cex=0.5)
##################################################################################################### 
####################### Step 6: TSNE to show subtype ############################################### 
##################################################################################################### 
load("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/Pancancer/TCGA-BRCA-BUR-PAN-TopVar.RData")
test<-input[,c(1,match("cg17780246",colnames(input)))]
boxplot(test$cg17780246~test$phen)