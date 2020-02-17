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
            cexRow = 0.5, cexCol = 0.5,
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
pdf("heatmap2.pdf")
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




############################################################################################################
### GSI and HEATMAP to DNA methylation beta matrix

setwd("/home/local/MFLDCLIN/guosa/hpc/project/thyroid/scRNA")

sam<-read.csv("https://raw.githubusercontent.com/Shicheng-Guo/thyroidscrna/master/extdata/scTCGA2019.csv")

# data$Manual_Type.according.to.SC.seq
# match(colnames(input),data$Patient_barcode)
# data$Patient_barcode[which(! data$Patient_barcode %in% colnames(input))]

load("/home/guosa/hpc/methylation/Pancancer/methdata.pancancer.RData")
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/GscTools.R")

id2phen4<-function(filename){
  library("stringr")
  as.array(str_extract(filename,"TCGA-[0-9|a-z|A-Z]*-[0-9|a-z|A-Z]*-[0-9]*"))
}

id2phen3<-function(filename){
  library("stringr")
  as.array(str_extract(filename,"TCGA-[0-9|a-z|A-Z]*-[0-9|a-z|A-Z]*"))
}

id2bin<-function(filename){
  library("stringr")
  filename<-as.array(str_extract(filename,"TCGA-[0-9|a-z|A-Z]*-[0-9|a-z|A-Z]*-[0-9]*"))
  as.numeric(lapply(strsplit(filename,"-"),function(x) x[4]))
}

input<-methdata
phen2<-id2bin(colnames(input))
input<-input[,which(phen2==1)]
phen3<-id2phen3(colnames(input))
input<-input[,match(sam$Patient_barcode,phen3)]
colnames(input)<-id2phen3(colnames(input))
newinput<-RawNARemove(input)
newinput<-RawZeroRemove(newinput)
colnames(newinput)<-colnames(input)

sctype<-sam[match(colnames(newinput),sam$Patient_barcode),]$Manual_Type.according.to.SC.seq
colnames(newinput)<-sctype
newinput<-newinput[,order(colnames(newinput))]
rlt<-GSIMethDecon(newinput)
xsel<-subset(rlt,GSI>0.43)
dim(xsel)
matrix<-newinput[match(xsel[order(xsel$group),1],rownames(newinput)),]

save(matrix,file="scRNA-thryoid-tcga-methylation.RData")

setwd("~/hpc/project/thyroid/scRNA")

library("gplots")
pdf("meth-heatmap2.pdf")
heatmap.2(scale(matrix),Rowv=NA,Colv=NA)
dev.off()

library("gplots")
pdf("meth-heatmap3.pdf")
heatmap.2(scale(matrix),Rowv=F,Colv=F,keysize=1, density.info="none",trace="none")
dev.off()

pdf("meth-heatmap4.pdf")
HeatMap(scale(na.omit(matrix)),Rowv=T,Colv=F)
dev.off()

pdf("meth-heatmap5.pdf")
HeatMap(matrix,Rowv=F,Colv=F,)
dev.off()

pdf("meth-heatmap6.pdf")
HeatMap(scale(matrix),Rowv=F,Colv=F)
dev.off()



