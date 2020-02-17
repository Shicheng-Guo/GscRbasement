##############################################################################################################
# The script was used to do co-expression analysis between miR32A and 104 TFs within miR32A promoter regions
# Shihcheng.Guo@Gmail.com
# 2020/01/22
# miR23A, miR27A, miR24-2
##############################################################################################################

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

RawNARemove<-function(data,missratio=0.3){
  threshold<-(missratio)*ncol(data)
  NaRaw<-which(apply(data,1,function(x) sum(is.na(x))>=threshold))
  zero<-which(apply(data,1,function(x) all(x==0))==T)
  NaRAW<-c(NaRaw,zero)
  if(length(NaRAW)>0){
    output<-data[-NaRAW,]
  }else{
    output<-data;
  }
  output
}

RawZeroRemove<-function(data,missratio=0.5){
  threshold<-(missratio)*ncol(data)
  NaRaw<-which(apply(data,1,function(x) sum(is.na(x))>=threshold))
  zero<-which(apply(data,1,function(x) sum(x==0)>=threshold))
  NaRAW<-c(NaRaw,zero)
  if(length(NaRAW)>0){
    output<-data[-NaRAW,]
  }else{
    output<-data;
  }
  output
}


load("~/hpc/methylation/Pancancer/RNA-seq/rnaseqdata.pancancer.RData")
phen1=read.table("~/hpc/methylation/TCGA-clinical-11093.tsv",header = T,sep="\t")
phen2=read.table("~/hpc/methylation/Pancancer/RNA-seq/File_metadata2.txt",header = T,sep="\t")

phen<-data.frame(phen2,phen1[match(phen2$cases.0.case_id,phen1$case_id),])
phen$file_name=gsub(".gz","",phen$file_name)
head(phen1)
head(phen2)
head(phen)

colnames(rnaseqdata)<-unlist(lapply(strsplit(colnames(rnaseqdata),"/"),function(x) x[2]))
phen<-phen[match(colnames(rnaseqdata),phen$file_name),]

phen1<-id2phen4(phen$cases.0.samples.0.submitter_id)
phen2<-id2phen3(phen$cases.0.samples.0.submitter_id)
phen3<-id2bin(phen$cases.0.samples.0.submitter_id)

phen$bin=phen3

include<-which(c(phen$bin==1 | phen$bin==11))
phen<-phen[include,]
input<-rnaseqdata[,include]
phen$id=id2phen4(phen$cases.0.samples.0.submitter_id)
dim(phen)
dim(input)
input[1:5,1:5]
phen[1:5,1:5]
colnames(input)<-phen$id
Seq<-paste(phen$project_id,phen$bin,sep="-")

load("~/hpc/project/TCGA/pancancer/miRNA/data/TCGA-Pancancer.miRNAseq.RData")

rnaseqdata<-input

miRNA[1:5,1:4]

ENSG2Symbol<-function(ENSG){
  db<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/ENSG.ENST.ENSP.Symbol.hg19.bed",sep="\t")
  ENSG<-unlist(lapply(strsplit(ENSG,split="[.]"),function(x) x[1]))
  Symbol<-db[match(as.character(ENSG),db$V8),4]
  return(Symbol)
}

symbol<-as.character(ENSG2Symbol(rownames(rnaseqdata)))
tfbs<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/tcga/master/extdata/gene/miR23A/miR23A-TF.txt",head=T,sep="\t",as.is=T)
tfbsrnaseq<-rnaseqdata[unique(na.omit(match(tfbs[,1],symbol))),]
tfbsmiRNAseq<-miRNA[,na.omit(match(colnames(tfbsrnaseq),colnames(miRNA)))]
tfbsrnaseq<-tfbsrnaseq[,match(colnames(tfbsmiRNAseq),colnames(tfbsrnaseq))]
rownames(tfbsrnaseq)=as.character(ENSG2Symbol(rownames(tfbsrnaseq)))
tfbsrnaseq[1:5,1:4]
tfbsmiRNAseq[1:5,1:4]

tfbsrnaseq<-log(tfbsrnaseq+1,2)
y<-log(tfbsmiRNAseq[match(c('hsa-mir-23a',"hsa-mir-27a","hsa-mir-24-2"),rownames(tfbsmiRNAseq)),]+1,2)
matrix<-data.frame(t(y),t(tfbsrnaseq),check.names = F)

fit<-c()
for(i in 1:nrow(tfbsrnaseq)){
  temp<-summary(lm(y~tfbsrnaseq[i,]))
  fit<-rbind(fit,temp$coefficients[2,])
}
rownames(fit)<-rownames(tfbsrnaseq)
fit<-fit[order(fit[,4]),]
write.csv(fit,file="miR32A-TFs.coexpression.csv",quote=F)

sel<-match('hsa-mir-32',rownames(tfbsmiRNAseq))
tfbsmiRNAseqmiR32A<-tfbsmiRNAseq[(sel-2):(sel+2),]
write.csv(tfbsmiRNAseqmiR32A,file="miRNA.csv",quote=F)
write.csv(tfbsrnaseq,file="mRNAseq.csv",quote=F)

x<-tfbsrnaseq[match("RUNX1T1",rownames(tfbsrnaseq)),]

pdf("RUNX1T1.pdf")
plot(y~x,pch=16,cex=0.3,col="red",xlab="RUNX1T1",ylab="miR23A")
dev.off()

x<-tfbsrnaseq[match("GTF2B",rownames(tfbsrnaseq)),]

pdf("GTF2B.pdf")
plot(y~x,pch=16,cex=0.3,col="red",xlab="GTF2B",ylab="miR23A")
dev.off()


##########################
# miR23A, miR27A, miR24-2
source("https://raw.githubusercontent.com/Shicheng-Guo/tcga/master/extdata/gene/miR23A/HeatMap.R")
x1<-rownames(miRNA)[grep("-23a",rownames(miRNA))]
x2<-rownames(miRNA)[grep("-27a",rownames(miRNA))]
x3<-rownames(miRNA)[grep("-24-2",rownames(miRNA))]
matrix<-t(miRNA[match(c(x1,x2,x3),rownames(miRNA)),])
cor(matrix)
pdf("hsa-miR23A-h3.pdf")
heatmap.3(cor(matrix))
dev.off()
pdf("hsa-miR23A-Heat.pdf")
HeatMap(cor(matrix))
dev.off()
tfbsrnaseq<-log(tfbsrnaseq+1,2)
y<-log(tfbsmiRNAseq[match(c('hsa-mir-32',"hsa-mir-27a","hsa-mir-24-2"),rownames(tfbsmiRNAseq)),]+1,2)
matrix<-data.frame(t(y),t(tfbsrnaseq),check.names = F)
pdf("hsa-miR-mRNA-FullHeatMap.pdf")
HeatMap(cor(matrix))
dev.off()
write.csv(cor(matrix),file="miR23A.correlation.csv",quote=F)

##########################
# miR23A, miR27A, miR24-2

load("~/hpc/methylation/Pancancer/RNA-seq/rnaseqdata.pancancer.RData")
phen1=read.table("~/hpc/methylation/TCGA-clinical-11093.tsv",header = T,sep="\t")
phen2=read.table("~/hpc/methylation/Pancancer/RNA-seq/File_metadata2.txt",header = T,sep="\t")
phen<-data.frame(phen2,phen1[match(phen2$cases.0.case_id,phen1$case_id),])
phen$file_name=gsub(".gz","",phen$file_name)
head(phen1)
head(phen2)
head(phen)
colnames(rnaseqdata)<-unlist(lapply(strsplit(colnames(rnaseqdata),"/"),function(x) x[2]))
phen<-phen[match(colnames(rnaseqdata),phen$file_name),]
phen1<-id2phen4(phen$cases.0.samples.0.submitter_id)
phen2<-id2phen3(phen$cases.0.samples.0.submitter_id)
phen3<-id2bin(phen$cases.0.samples.0.submitter_id)
phen$bin=phen3
include<-which(c(phen$bin==1 | phen$bin==11))
phen<-phen[include,]
input<-rnaseqdata[,include]
phen$id=id2phen4(phen$cases.0.samples.0.submitter_id)
dim(phen)
dim(input)
input[1:5,1:5]
phen[1:5,1:5]
colnames(input)<-phen$id
Seq<-paste(phen$project_id,phen$bin,sep="-")
input[1:5,1:5]
rownames(input)=as.character(ENSG2Symbol(rownames(input)))
input[1:5,1:5]
input<-t(log(input+1,2))
tfbs<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/ferroptosis/master/NCOA4.tfbs.txt",sep="\t")
matrix<-input[,match(tfbs[,1],colnames(input))]
pdf("NCOA4-mRNA-FullHeatMap.V3.pdf")
HeatMap(cor(matrix),cexRow=0.2, cexCol=0.2)
dev.off()
write.csv(cor(matrix),file="NCOA4.correlation.v2.csv",quote=F)

matrix<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/ferroptosis/master/NCOA4.TFBS.coexpressionMatrix.txt")






matrix<-read.csv("https://raw.githubusercontent.com/Shicheng-Guo/ferroptosis/master/extdata/NCOA4-blocktcga.pancancer.meta.dge.smd.pvalue.csv")
matrix<-matrix[-3,]
matrix<-head(matrix,n=30)
matrix$pval.random=-log(matrix$pval.random,10)
matrix<-matrix[order(matrix$pval.random,decreasing = F),]
head(matrix)
par(las=1,cex=0.75)
barplot(matrix$pval.random,col=2:nrow(matrix),horiz = T,names.arg =matrix$symbol)

x1<-unlist(lapply(strsplit(as.character(sema[,1]),"_"),function(x) x[2]))
x2<-unlist(lapply(strsplit(as.character(broad[,1]),"-"),function(x) x[1]))

sum(x1 %in% x2)



polyphen<-read.table("polyphen.txt",head=T,sep="\t")
sift<-read.table("sift.txt",head=T,sep="\t")
head(polyphen)
head(sift)
out<-merge(polyphen,sift,by="Variation.ID")
head(out)
rlt<-data.frame(out$Variation.ID,out$Prediction.x,out$Prediction.y,out$Wild.AA.x,out$AA.Position.x,out$Mutant.AA.x)
output<-unique(subset(rlt,out.Prediction.y=="Deleterious" & out.Prediction.x=="Probably Damaging"))
head(output)
write.csv(unique(output),file="RALGPS1.csv",quote=F,row.names=F)

data<-read.table("avinput.txt.hg19_multianno.txt",head=T,sep="\t",as.is=T)

RiskAAC<-function(data){
  risk1<-grep("D",data$SIFT_pred)
  risk2<-grep("D|P",data$Polyphen2_HDIV_pred)
  risk3<-grep("D|P",data$Polyphen2_HVAR_pred)
  risk4<-grep("D",data$LRT_pred)
  risk5<-grep("D",data$MutationTaster_pred)
  risk6<-grep("H|M",data$MutationAssessor_pred)
  risk7<-grep("D",data$FATHMM_pred)
  risk8<-grep("D",data$PROVEAN_pred)
  risk9<-grep("D",data$MetaSVM_pred)
  risk10<-grep("D",data$MetaLR_pred)
  risk11<-grep("D",data$fathmm.MKL_coding_pred)
  risk12<-grep("D",data$M.CAP_pred)
  rlt=c(risk1,risk2,risk3,risk4,risk5,risk6,risk7,risk8,risk9,risk10,risk11,risk12)
  return(rlt)
}

data$nsyn<-table(RiskAAC(data))
write.table(data,file="annovar.output.txt",quote=F,sep="\t",col.names = T,row.names = F)

