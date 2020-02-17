setwd("~/hpc/methylation/Pancancer/RNA-seq")
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/GscTools.R")
library("metafor")
library("meta")
library("metacor")

load("~/hpc/methylation/Pancancer/RNA-seq/rnaseqdata.pancancer.RData")


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

GSIMethDecon<-function(data){
  data<-data.matrix(data)
  group=names(table(colnames(data)))
  index=colnames(data)
  gsi<-gmaxgroup<-avebase<-c()
  for(i in 1:nrow(data)){
  gsit<-0
  gmax<-names(which.max(tapply(as.numeric(data[i,]),index,function(x) mean(x,na.rm=T))))
  for(j in 1:length(group)){
    tmp<-(1-10^(mean(data[i,][which(index==group[j])],na.rm=T))/10^(mean(data[i,][which(index==gmax)],,na.rm=T)))/(length(group)-1)
    gsit<-gsit+tmp
  }
  ave<-tapply(data[i,], index, function(x) mean(x,na.rm=T))
  gmaxgroup<-c(gmaxgroup,gmax)
  gsi<-c(gsi,gsit)
  avebase<-rbind(avebase,ave)
  }
  rlt=data.frame(region=rownames(data),group=gmaxgroup,GSI=gsi,AVE=avebase)
  return(rlt)
}

gsi<-function(data){
  group=names(table(colnames(data)))
  index=colnames(data)
  GSI<-c()
  gmaxgroup<-c()
  for(i in 1:nrow(data)){
    gsit<-0
    gmax<-names(which.max(tapply(as.numeric(data[i,]),index,function(x) mean(x,na.rm=T))))
    if(length(gmax)<1){print(data[i,])}
    for(j in 1:length(group)){
      tmp<-(1-10^(mean(na.omit(as.numeric(data[i,which(index==group[j])])),na.rm=T))/10^(mean(na.omit(as.numeric(data[i,which(index==gmax)])))))/(length(group)-1)
      gsit<-c(gsit,tmp)
    }
    gmaxgroup<-c(gmaxgroup,gmax)
    GSI<-c(GSI,sum(gsit,na.rm=T))  # debug for GSI=NAN
  }
  rlt=data.frame(region=rownames(data),group=gmaxgroup,GSI=GSI)
  return(rlt)
}

HeatMap<-function(data,Rowv=T,Colv=T){
  
#  note: this function include correlation based heatmap (pearson or spearman)
#  data: row is gene and column is sample
#  colname and rowname cannot be NULL  
#  Usage example:
#  test<- matrix(runif(100),nrow=20)
#  colnames(test)=c("A","A","A","B","B")
#  rownames(test)=paste("Gene",1:20,sep="")
#  HeatMap(test)
  
  library("gplots")
  colors <- colorpanel(75,"midnightblue","mediumseagreen","yellow") 
  colors <-bluered(75)
  
  sidecol<-function(x){
    x<-as.numeric(as.factor(x))
    col<-rainbow(length(table(colnames(data))))
    sapply(x,function(x) col[x])
  }
  
  Hclust=function(x){hclust(x,method="complete")}
  Distfun=function(x){as.dist((1 - cor(t(x),method = "spearman")))}
  
  ColSideColors=sidecol(colnames(data))
  
  heatmap.2(data,trace="none", 
            hclust=Hclust,
            distfun=Distfun, 
            cexRow = 1, cexCol = 1,
            ColSideColors=ColSideColors,
            density.info="none",col=colors,
            Colv=Colv,Rowv = Rowv,
            keysize=0.9, margins = c(5, 10)
            )
}


TCGAProjects=c("BLCA","BRCA","CESC","CHOL","COAD","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC","LUAD","LUSC","PAAD","PCPG","PRAD","READ","SARC","STAD","THCA","THYM","UCEC")
TCGAProjects=c("THCA")
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

idx<-which(phen$phen2==1)
phen<-phen[idx,]
input<-rnaseqdata[,idx]

idx<-which(phen$pid %in% paste("TCGA-",TCGAProjects,sep=""))
phen<-phen[idx,]
input<-input[,idx]
input[1:5,1:5]

input<-log(input+1,2)
input<-RawNARemove(input)
input<-RawZeroRemove(input)
input[1:5,1:5]

colnames(input)<-phen$cases.0.submitter_id
input[1:5,1:5]

############################################################################################################
### GSI and HEATMAP (scRNA-seq)
setwd("/home/local/MFLDCLIN/guosa/hpc/project/thyroid/scRNA")
data<-read.csv("https://raw.githubusercontent.com/Shicheng-Guo/thyroidscrna/master/extdata/scTCGA2019.csv")
# data$Manual_Type.according.to.SC.seq
# match(colnames(input),data$Patient_barcode)
# data$Patient_barcode[which(! data$Patient_barcode %in% colnames(input))]
newinput<-input[,which(colnames(input) %in% data$Patient_barcode)]
sctype<-data[match(colnames(newinput),data$Patient_barcode),]$Manual_Type.according.to.SC.seq
colnames(newinput)<-sctype
newinput<-newinput[,order(colnames(newinput))]
rlt<-GSIMethDecon(newinput)
xsel<-subset(rlt,GSI>0.99)
matrix<-newinput[match(xsel[order(xsel$group),1],rownames(newinput)),]

library("gplots")
pdf("heatmap2.pdf")
heatmap.2(scale(matrix),Rowv=NA,Colv=NA)
dev.off()

library("gplots")
pdf("heatmap3.pdf")
heatmap.2(scale(matrix),Rowv=F,Colv=F,keysize=1, density.info="none",trace="none")
dev.off()

pdf("heatmap4.pdf")
HeatMap(scale(matrix),Rowv=T,Colv=F)
dev.off()

pdf("heatmap6.pdf")
HeatMap(matrix,Rowv=F,Colv=F)
dev.off()

pdf("heatmap7.pdf")
HeatMap(scale(matrix),Rowv=F,Colv=F)
dev.off()





xxxv<-read.csv("https://raw.githubusercontent.com/Shicheng-Guo/ferroptosis/master/ferroptosis.genelist.csv",head=F)
xxxv<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/ferroptosis/master/tsg.positivecontrol.txt",head=F)
xxxv<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/esophageal/master/phase1.genelist.txt",head=F)

ENSG<-Symbol2ENSG(as.character(xxxv[,1]))
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
head(rlt[order(rlt$pval),])
write.table(rlt,file="/home/local/MFLDCLIN/guosa/hpc/project/ferroptosis/ferroptosis.meta.table.pvalue.txt",sep="\t",quote=F,col.names = NA,row.names = T)

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
				 
				 












