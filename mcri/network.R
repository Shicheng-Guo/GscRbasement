##############################################################################################
##############################################################################################
##############################################################################################
manifest="gdc_manifest.2019-03-09.txt"
x=read.table(manifest,header = T)
manifest_length= nrow(x)
id= toString(sprintf('"%s"', x$id))
Part1= '{"filters":{"op":"in","content":{"field":"files.file_id","value":[ '
Part2= '] }},"format":"TSV","fields":"file_id,file_name,cases.submitter_id,cases.case_id,data_category,data_type,cases.samples.tumor_descriptor,cases.samples.tissue_type,cases.samples.sample_type,cases.samples.submitter_id,cases.samples.sample_id,cases.samples.portions.analytes.aliquots.aliquot_id,cases.samples.portions.analytes.aliquots.submitter_id","size":'
Part3= paste0("\"",manifest_length, "\"", "}")
Sentence= paste(Part1,id,Part2,Part3, collapse=" ")
write.table(Sentence,"Payload.txt",quote=F,col.names=F,row.names=F)
system("curl --request POST --header \"Content-Type: application/json\" --data @Payload.txt \"https://api.gdc.cancer.gov/files\" > file_metadata.txt")

##############################################################################################
##############################################################################################
##############################################################################################
## 

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

setwd("~/hpc/methylation/Pancancer/RNA-seq")

load("rnaseqdata.pancancer.RData")
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

include<-which(c(phen$bin==1 |phen$bin==11))
phen<-phen[include,]
input<-rnaseqdata[,include]
phen$id=id2phen4(phen$cases.0.samples.0.submitter_id)
dim(phen)
dim(input)
head(phen)
colnames(input)<-phen$id
Seq<-paste(phen$project_id,phen$bin,sep="-")
pid<-phen$project_id
input[1:5,1:5]
phen[1:5,1:5]
dim(input)
dim(phen)
##############################################################################################
##############################################################################################
##############################################################################################
library("metafor")

data<-input
i=500
Seq<-paste(phen$project_id,phen$bin,sep="-")
mean<-tapply(as.numeric(data[i,]),Seq,function(x) mean(x,na.rm=T))
sd<-tapply(as.numeric(data[i,]),Seq,function(x) sd(x,na.rm=T))
num<-tapply(as.numeric(data[i,]),Seq,function(x) length(x))

exclude<-names(which(table(unlist(lapply(strsplit(names(mean),"-"),function(x) x[2])))<2))
exclude <-grep(paste(exclude,collapse="|"),phen$project_id)

phen<-phen[-exclude,]
input<-input[,-exclude]
phen$id=id2phen4(phen$cases.0.samples.0.submitter_id)
colnames(input)<-phen$id
dim(phen)
dim(input)
input[1:5,1:5]
phen[1:5,1:5]

input<-log(input+1,2)
input<-RawNARemove(input)
input<-RawZeroRemove(input)
Seq<-paste(phen$project_id,phen$bin,sep="-")

rlt<-c()
coll<-c()
for(i in 1:nrow(input)){
  print(i)
  mean<-tapply(as.numeric(input[i,]),Seq,function(x) mean(x,na.rm=T))
  sd<-tapply(as.numeric(input[i,]),Seq,function(x) sd(x,na.rm=T))
  num<-tapply(as.numeric(input[i,]),Seq,function(x) length(x))
  m1i=mean[seq(1,length(mean),by=2)]
  m2i=mean[seq(2,length(mean),by=2)]
  sd1i=sd[seq(1,length(mean),by=2)]
  sd2i=sd[seq(2,length(mean),by=2)]
  n1i=num[seq(1,length(mean),by=2)]
  n2i=num[seq(2,length(mean),by=2)]
  Source<-unlist(lapply(strsplit(names(m1i),"_"),function(x) x[1]))
  output<-data.frame(cbind(n1i,m1i,sd1i,n2i,m2i,sd2i))
  output$source=Source
  output<-na.omit(output)
  es<-escalc(m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i,measure="MD",data=output)
  md <- rma(es,slab=source,method = "REML", measure = "SMD",data=output)
  rlt<-rbind(rlt,c(i,md$beta,md$pval,md$ci.lb,md$ci.ub,md$I2,md$tau2))
  coll<-c(coll,i)
}
rownames(rlt)<-rownames(input)[coll]
colnames(rlt)<-c("idx","beta","pval","cilb","ciub","i2","tau2")
rlt<-data.frame(rlt)
head(rlt[order(rlt$pval),])

i=8761
print(i)
mean<-tapply(as.numeric(input[i,]),Seq,function(x) mean(x,na.rm=T))
sd<-tapply(as.numeric(input[i,]),Seq,function(x) sd(x,na.rm=T))
num<-tapply(as.numeric(input[i,]),Seq,function(x) length(x))
m1i=mean[seq(1,length(mean),by=2)]
m2i=mean[seq(2,length(mean),by=2)]
sd1i=sd[seq(1,length(mean),by=2)]
sd2i=sd[seq(2,length(mean),by=2)]
n1i=num[seq(1,length(mean),by=2)]
n2i=num[seq(2,length(mean),by=2)]
Source<-unlist(lapply(strsplit(names(m1i),"_"),function(x) x[1]))
output<-data.frame(cbind(n1i,m1i,sd1i,n2i,m2i,sd2i))
output$source=Source
output<-na.omit(output)
es<-escalc(m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i,measure="MD",data=output)
md <- rma(es,slab=source,method = "REML", measure = "SMD",data=output)
print(md)
dev.off()


install.packages('biomaRt')
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- df$genes
df<-df[,-4]
G_list <- getBM(filters= "ensembl_peptide_id", attributes= c("ensembl_peptide_id","hgnc_symbol"),values=genes,mart= mart)
merge(df,G_list,by.x="gene",by.y="ensembl_peptide_id")

##############################################################################################
############### Co-expression Network Analysis (BIRC10) ################################
##############################################################################################
library("CEMiTool")

setwd("/mnt/bigdata/Genetic/Projects/shg047/methylation/Pancancer/RNA-seq")

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


ENSG<-Symbol2ENSG(c("TP53","KRAS","TTN","MUC16","CSMD3","PIK3CA","CDK2"))
xgene<-c("ENSG00000213754",as.character(ENSG[,2]))
xsel<-unlist(lapply(xgene,function(x) grep(x,rownames(input))))
name<-as.character(ENSG2Symbol(xgene))
name<-c("PANC754","TP53","KRAS","TTN","MUC16","CSMD3","PIK3CA","CDK2")
temp<-input[xsel,]
rownames(temp)<-name
value<-cor(t(temp))
value
pdf("../../PANC754.TP53.heatmap.pdf")
HeatMap(cor(t(temp)),cexRow=0.75,cexCol=0.75)
dev.off()
write.table(value,file="../../PANC754.TP53.correlation.txt",sep="\t",quote=F,col.names = NA,row.names = T)
P<-c()
corvalue<-c()
for(i in 2:nrow(temp)){
  xx<-cor.test(temp[1,],temp[i,])
  corvalue<-c(corvalue,xx$estimate)
  P<-c(P,xx$p.value)
}

rlt<-data.frame(corvalue,P)
rownames(rlt)<-rownames(temp)[2:nrow(temp)]
write.table(rlt,file="../../PANC754.TP53.co-expression.pvalue.txt",sep="\t",col.names = NA,row.names = T,quote=F)



############################################
xgene<-c("ENSG00000213754")
xsel<-unlist(lapply(xgene,function(x) grep(x,rownames(input))))
xsel
P<-c()
corvalue<-c()
for(i in 1:nrow(input)){
  xx<-cor.test(input[xsel,],input[i,])
  corvalue<-c(corvalue,xx$estimate)
  P<-c(P,xx$p.value)
}

rlt<-data.frame(corvalue,P)
rownames(rlt)<-rownames(input)
name<-ENSG2Symbol(rownames(rlt))
rlt<-na.omit(data.frame(name,rlt))
out<-na.omit(rlt)
padj<-p.adjust(out$P,method="fdr")
out<-data.frame(out,padj)
padj.bonferroni<-p.adjust(out$P,method="bonferroni")
out<-data.frame(out,padj.bonferroni)
sum(out$padj.bonferroni<0.05)

db<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/ENSG.ENST.ENSP.Symbol.hg19.bed",sep="\t")

output<-data.frame(db[match(out$name,db$V4),],out)

setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/Pancancer/RNA-seq")
load("output.RData")
chr<-gsub("chr","",output$V1)
chr[chr=="X"]<-23
chr[chr=="Y"]<-24
chr<-as.numeric(chr)
cm1<-na.omit(data.frame(SNP=output$name,Chromosome=chr,Position=output$V2,P=output$P,padj=output$padj))
cm2<-na.omit(data.frame(SNP=output$name,chr=chr,pos=output$V2,R=output$corvalue,P=output$P,padj=output$padj))

library("CMplot")
CMplot(cm, plot.type="m", band=0, LOG10=FALSE, ylab="Correlation Coefficient",ylim=c(-0.2,0.6),
       threshold.lty=2, threshold.lwd=1, threshold.col="red", amplify=FALSE,cex=0.6,
       chr.den.col=NULL, file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

CMplot(cm1,plot.type="m",LOG10=TRUE,threshold=NULL,file="jpg",memo="",dpi=300,ylim=c(0,350),
       file.output=TRUE,verbose=TRUE,width=14,height=6)

dim(subset(cm,R>0.3))
dim(subset(cm,R< -0.12))
write.table(subset(cm,R>0.3 | R< -0.12),file="PANC754.signficant.correlated.genes.txt",sep="\t",col.names = NA,row.names = T,quote=F)



sym<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/PANC/master/extdata/PANC754/alignment.txt",head=F)
xint<-Symbol2ENSG(unique(as.character(sym[,1])))
xgene<-c("ENSG00000213754.2",as.character(xint[,2]))
xsel<-unlist(lapply(xgene,function(x) grep(x,rownames(input))))

temp<-input[xsel,]
temp<-log(temp+1,2)

name<-as.character(ENSG2Symbol(rownames(temp)))
name[1]<-"PANC754"
rownames(temp)<-name
value<-cor(t(temp))
value

HeatMap<-function(data,Rowv=T,Colv=T,cexRow=0.5,cexCol=0.5){
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
            cexRow = cexRow, cexCol = cexCol,
            ColSideColors=ColSideColors,
            density.info="none",col=colors,
            Colv=Colv,Rowv = Rowv,
            keysize=0.9, margins = c(5, 10)
  )
} 

library("gplots")
pdf("../../heatmap.pdf")
HeatMap(cor(t(temp)),cexRow=0.75,cexCol=0.75)
dev.off()
write.table(value,file="../../PANC754.correlation.txt",sep="\t",quote=F,col.names = NA,row.names = T)

P<-c()
corvalue<-c()
for(i in 2:nrow(temp)){
  xx<-cor.test(temp[1,],temp[i,])
  corvalue<-c(corvalue,xx$estimate)
  P<-c(P,xx$p.value)
}

rlt<-data.frame(corvalue,P)
rownames(rlt)<-rownames(temp)[2:nrow(temp)]
write.table(rlt,file="../../co-expression.pvalue.txt",sep="\t",col.names = NA,row.names = T,quote=F)

out<-c()
name<-c()
for(i in unique(pid)){
  temp<-input[,which(pid==i)]
  print(c(i,dim(temp)))
  if(ncol(data.frame(temp))>2){
  temp<-log(temp+1,2)
  temp<-temp[unlist(lapply(xgene,function(x) grep(x,rownames(temp)))),]
  fit<-cor.test(x=t(temp)[,1],y=t(temp)[,2])
  out<-rbind(out,c(fit$estimate,fit$p.value))
  name<-c(name,i)
  }
}
rownames(out)<-name
colnames(out)<-c("R","Pvalues")
write.csv(out,file="correlation.PANC754-KPNA4.cancer.txt")



xgene<-c("ENSG00000104643.8","ENSG00000213754.2")
out<-c()
name<-c()
for(i in unique(pid)){
  temp<-input[,which(pid==i)]
  print(c(i,dim(temp)))
  if(ncol(data.frame(temp))>2){
    temp<-log(temp+1,2)
    temp<-temp[unlist(lapply(xgene,function(x) grep(x,rownames(temp)))),]
    fit<-cor.test(x=t(temp)[,1],y=t(temp)[,2])
    out<-rbind(out,c(fit$estimate,fit$p.value))
    name<-c(name,i)
  }
}
rownames(out)<-name
colnames(out)<-c("R","Pvalues")
write.csv(out,file="correlation.PANC754-MTMR9.cancer.txt")




args = commandArgs(trailingOnly=TRUE)
meta.file=as.character(args[1])
output.file=as.character(args[2])
input<-read.table(meta.file,head=T,sep="\t")
table<-data.frame(table,OR=exp(table$bcac_onco_icogs_gwas_beta))
db<-read.table("~/hpc/db/hg19/m6ASNP.db147.hg19.uni.txt",head=T)
head(db)
newdb<-db[db$Rs_ID %in% table$phase3_1kg_id,]
out<-data.frame(newdb,table[match(newdb$Rs_ID,table$phase3_1kg_id),])
head(out)
write.table(out,file=output.file,sep="\t",quote=F,row.names =F,col.names = T)





