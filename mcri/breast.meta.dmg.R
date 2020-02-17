source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/GscTools.R")
library("meta")
library("metafor")
library("survival")
library("survminer")

Symbol2ENSG<-function(Symbol){
  db<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/ENSG.ENST.ENSP.Symbol.hg19.bed",sep="\t")
  ENSG<-as.character(db[match(Symbol,db$V4),8])
  ENSG<-na.omit(data.frame(Symbol,ENSG))
  return(ENSG)
}
ENSG2Symbol<-function(ENSG){
  db<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/ENSG.ENST.ENSP.Symbol.hg19.bed",sep="\t")
  ENSG<-unlist(lapply(strsplit(ENSG,split="[.]"),function(x) x[1]))
  Symbol<-db[match(as.character(ENSG),db$V8),4]
  return(Symbol)
}

ensg2bed<-function(ENSG){
  db<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/hg19/ENSG.ENST.hg19.txt",as.is=T,head=F)
  ENSG<-unlist(lapply(strsplit(ENSG,split="[.]"),function(x) x[1]))
  bed<-unique(db[db$V5 %in% as.character(ENSG),c(1,2,3,5)])
  return(bed)
}


source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/GscTools.R")
library("metafor")
library("meta")
library("metacor")

load("~/hpc/methylation/Pancancer/RNA-seq/rnaseqdata.pancancer.RData")


TCGAProjects=c("BLCA","BRCA","CESC","CHOL","COAD","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC","LUAD","LUSC","PAAD","PCPG","PRAD","READ","SARC","STAD","THCA","THYM","UCEC")
phen1=read.table("~/hpc/methylation/TCGA-clinical-11093.tsv",header = T,sep="\t")
phen2=read.table("~/hpc/methylation/Pancancer/RNA-seq/File_metadata2.txt",header = T,sep="\t")
head(phen1)
head(phen2)
colnames(rnaseqdata)<-unlist(lapply(strsplit(colnames(rnaseqdata),"/"),function(x) x[2]))
phen<-data.frame(phen2,phen1[match(phen2$cases.0.case_id,phen1$case_id),])
phen$file_name=gsub(".gz","",phen$file_name)

# prepare phenotype information
phen<-phen[match(colnames(rnaseqdata),phen$file_name),]
phen$phen4<-id2phen4(phen$cases.0.samples.0.submitter_id)
phen$phen3<-id2phen3(phen$cases.0.samples.0.submitter_id)
phen$phen2<-id2bin(phen$cases.0.samples.0.submitter_id)
phen$pid<-phen$project_id
head(phen)

idx<-which(phen$phen2==1 | phen$phen2==11)
phen<-phen[idx,]
input<-rnaseqdata[,idx]

idx<-which(phen$pid %in% paste("TCGA-",TCGAProjects,sep=""))
phen<-phen[idx,]
input<-input[,idx]
input[1:5,1:5]

input<-log(input+1,2)
input<-RawNARemove(input)
input<-RawZeroRemove(input)

xxxv<-read.csv("https://raw.githubusercontent.com/Shicheng-Guo/ferroptosis/master/ferroptosis.genelist.csv",head=F)
xxxv<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/ferroptosis/master/tsg.positivecontrol.txt",head=F)
xxxv<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/esophageal/master/phase1.genelist.txt",head=F)
xxxv<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/breast/master/target.txt",head=T)
xxxv<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/cholangiocarcinoma/master/cholangiocarcinoma.hg19.bed",head=T,as.is=T,sep="\t")

setwd("~/hpc/methylation/Pancancer/RNA-seq")


ENSG<-Symbol2ENSG(as.character(xxxv[,4]))
xgene<-c(as.character(ENSG[,2]))
ii<-unlist(lapply(xgene,function(x) grep(x,rownames(input))))

Seq<-paste(phen$project_id,phen$phen2,sep="-")
rlt<-c()
coll<-c()

for(i in ii){
  mean<-tapply(as.numeric(input[i,]),Seq,function(x) mean(x,na.rm=T))
  sd<-tapply(as.numeric(input[i,]),Seq,function(x) sd(x,na.rm=T))
  num<-tapply(as.numeric(input[i,]),Seq,function(x) length(x))
  m1i=mean[seq(1,length(mean),by=2)]
  m2i=mean[seq(2,length(mean),by=2)]
  sd1i=sd[seq(1,length(mean),by=2)]
  sd2i=sd[seq(2,length(mean),by=2)]
  n1i=num[seq(1,length(mean),by=2)]
  n2i=num[seq(2,length(mean),by=2)]
  Source<-unlist(lapply(strsplit(names(m1i),"-"),function(x) x[2]))
  output<-data.frame(cbind(n1i,m1i,sd1i,n2i,m2i,sd2i))
  output$source=Source
  output<-na.omit(output)
  es<-escalc(m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i,measure="MD",data=output)
  md <- rma(es,slab=source,method = "REML", measure = "SMD",data=output)
  rlt<-rbind(rlt,c(i,md$beta,md$pval,md$ci.lb,md$ci.ub,md$I2,md$tau2))
  coll<-c(coll,i)
  m<-metagen(yi,seTE=vi,data = es,
             comb.fixed = TRUE,
             comb.random = TRUE,
             prediction=F,
             sm="SMD")

  Symbol<-ENSG2Symbol(rownames(input)[i])
  print(c(i,as.character(Symbol)))
  pdf(paste("/home/local/MFLDCLIN/guosa/hpc/project/ferroptosis/",Symbol,"-",rownames(input)[i],".SMD.PANC.pdf",sep=""))
  forest(m,leftlabs = Source,
         lab.e = "Intervention",
         pooled.totals = FALSE,
         smlab = "",studlab=Source,
         text.random = "Overall effect",
         print.tau2 = FALSE,
         col.diamond = "blue",
         col.diamond.lines = "black",
         col.predict = "red",
         print.I2.ci = TRUE,
         digits.sd = 2,fontsize=8)
  dev.off()
}
rownames(rlt)<-ENSG2Symbol(rownames(input)[coll])
colnames(rlt)<-c("idx","beta","pval","cilb","ciub","i2","tau2")
rlt<-data.frame(rlt)
rlt<-rlt[order(rlt$pval),]
write.table(rlt,file="/home/local/MFLDCLIN/guosa/hpc/project/ferroptosis/breast.meta.table.pvalue.txt",sep="\t",quote=F,col.names = NA,row.names = T)
write.csv(rlt,file="/home/local/MFLDCLIN/guosa/hpc/project/ferroptosis/breast.meta.table.pvalue.csv",quote=F)


rownames(rlt)<-rownames(input)[coll]
bed<-ensg2bed(as.character(rownames(input)[coll]))
xsel<-unlist(apply(bed,1,function(x) grep(x[4],rownames(rlt))))
output<-data.frame(bed,rlt[xsel,])
library("CMplot")
cminput<-data.frame(SNP=output$V5,Chromosome=output$V1,Position=output$V2,trait1=output[,7])
CMplot(cminput,plot.type="b",ylim=20,LOG10=TRUE,threshold=NULL,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6)

write.table(cminput,file="ferroptosis.pancancer.meta.dge.txt",sep="\t",quote=F,row.name=T,col.names=NA)

# co-expression network
netgene<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/ferroptosis/master/NCOA4.tfbs.txt",head=F)
ENSG<-Symbol2ENSG(as.character(c("NCOA4",as.character(netgene[,1]))))
xgene<-as.character(ENSG[,2])

Rmatrix<-c()
for(j in 2:length(xgene)){
out<-c()
name<-c()
	for(i in unique(pid)){
	  temp<-input[,which(pid==i)]
	  print(c(i,dim(temp)))
	  if(ncol(data.frame(temp))>2){
	  temp<-log(temp+1,2)
	  temp<-temp[unlist(lapply(xgene[c(1,j)],function(x) grep(x,rownames(temp)))),]
	  fit<-cor.test(x=t(temp)[,1],y=t(temp)[,2])
	  out<-rbind(out,c(fit$estimate,fit$p.value))
	  name<-c(name,i)
	  }
	rownames(out)<-name
	colnames(out)<-c("R","Pvalues")
	}
	xxgene<-ENSG2Symbol(rownames(temp))
	filename<-paste("../../",xxgene[1],"-",xxgene[2],".tcga.pancancer.coexpression.csv",sep="")
	print(xxgene[2])
	write.csv(out,file=filename)
}


m.cor <- metacor(cor, n, data = cordata, studlab = cordata$Author, sm = "ZCOR", method.tau = "SJ")
				 
				 












