

data<-read.table("AHRR.hap.txt",row.names=1)
sumhap<-unique(as.factor(c(as.character(data[,1]),as.character(data[,2]))))

data<-rbind(data,data.frame(V2=sumhap,V3=sumhap))

cor2bed<-function(cor){
  cor<-as.character(cor)
  a<-unlist(lapply(strsplit(cor,split=c(":")),function(x) strsplit(x,"-")))
  bed<-matrix(a,ncol=3,byrow=T)
  return(data.frame(bed))
}


table(data)

table(data,dnn=sumhap)
rw=unlist(lapply(strsplit(rownames(table(data)),split=""),function(x) sum(x %in% "C")))
cw=unlist(lapply(strsplit(colnames(table(data)),split=""),function(x) sum(x %in% "C")))
rrw=tapply(rowSums(table(data)),rw,sum)
ccw=tapply(colSums(table(data)),cw,sum)
vv=rbind(data.frame(v=rrw,w=names(rrw)),data.frame(v=ccw,w=names(ccw)))
tapply(vv$v,vv$w,sum)
table(data)
write.table(table(data),file="AHRR.dip.txt",sep="\t",quote=F,row.names=T,col.names=NA)



data=read.table("Nash.mf.matrix.txt",head=T,row.names=1,check.names=F)
colnames(data)=unlist(lapply(strsplit(colnames(data),split="_"),function(x) paste(x[1],x[2],sep="_")))
colnames(data)
threshold=30
pdata1=data[,c(2:4,18)]
pdata1[,4]=pdata1[,4]*100
head(pdata1)
F1=which(rowMeans(pdata1)>threshold)
length(F1)
rownames(F1)
sum((data[F1,c(31)])>threshold,na.rm=T)/length(na.omit((data[F1,c(31)])>threshold))
sum((data[F1,c(33)])>threshold,na.rm=T)/length(na.omit((data[F1,c(33)])>threshold))
sum((data[F1,c(35)])>threshold,na.rm=T)/length(na.omit((data[F1,c(35)])>threshold))
sum((data[F1,c(1)])>threshold,na.rm=T)/length(na.omit((data[F1,c(1)])>threshold))
sum(rowMeans(data[F1,c(31,33,35)])>threshold,na.rm=T)/length(na.omit(rowMeans(data[F1,c(31,33,35)])>threshold))
pdata2=data[F1,30:35]
newdata=pdata2[which(pdata2[,6]<pdata2[,4] & pdata2[,4]<pdata2[,2] & pdata2[,2]>pdata2[,1] & pdata2[,4]>pdata2[,3] & pdata2[,6]>pdata2[,5]),]
colnames(newdata)=unlist(lapply(strsplit(colnames(newdata),split="_"),function(x) x[2]))
dim(newdata)
target=subset(newdata,newdata[,2]>60)
input=data[match(rownames(target),rownames(data)),]
colnames(input)=unlist(lapply(strsplit(colnames(input),split="_"),function(x) paste(x[1],x[2],sep="_")))

# LINE-1
data=read.table("LINE-1.mf.txt",head=T,row.names=1,check.names=F)
LINE=colMeans(data,na.rm=T)
names(LINE)=unlist(lapply(strsplit(colnames(data),split="_"),function(x) paste(x[1],x[2],sep="_")))
t.test(LINE[2:4],LINE[c(31,33,35)])

colnames(data)
threshold=30
pdata1=data[,c(2:4,18)]
pdata1[,4]=pdata1[,4]*100
head(pdata1)
F1=which(rowMeans(pdata1)>threshold)
length(F1)
rownames(F1)
sum((data[F1,c(31)])>threshold,na.rm=T)/length(na.omit((data[F1,c(31)])>threshold))
sum((data[F1,c(33)])>threshold,na.rm=T)/length(na.omit((data[F1,c(33)])>threshold))
sum((data[F1,c(35)])>threshold,na.rm=T)/length(na.omit((data[F1,c(35)])>threshold))
sum((data[F1,c(1)])>threshold,na.rm=T)/length(na.omit((data[F1,c(1)])>threshold))


# Plan B

data=read.table("Nash.mf.matrix.txt",head=T,row.names=1,check.names=F)
colnames(data)=unlist(lapply(strsplit(colnames(data),split="_"),function(x) paste(x[1],x[2],sep="_")))
colnames(data)
threshold=60
pdata1=data[,c(2:4,18)]
pdata1[,4]=pdata1[,4]*100
head(pdata1)
F1=which(rowMeans(pdata1)>threshold)
length(F1)
rownames(F1)
sum((data[F1,c(31)])>threshold,na.rm=T)/length(na.omit((data[F1,c(31)])>threshold))
sum((data[F1,c(33)])>threshold,na.rm=T)/length(na.omit((data[F1,c(33)])>threshold))
sum((data[F1,c(35)])>threshold,na.rm=T)/length(na.omit((data[F1,c(35)])>threshold))
sum((data[F1,c(1)])>threshold,na.rm=T)/length(na.omit((data[F1,c(1)])>threshold))
sum(rowMeans(data[F1,c(31,33,35)])>threshold,na.rm=T)/length(na.omit(rowMeans(data[F1,c(31,33,35)])>threshold))

pdata2=data[F1,30:35]
newdata=pdata2[which(pdata2[,2]>pdata2[,1] & pdata2[,4]>pdata2[,3] & pdata2[,6]>pdata2[,5]),]
colnames(newdata)=unlist(lapply(strsplit(colnames(newdata),split="_"),function(x) x[2]))
dim(newdata)
newdata=subset(newdata,newdata[,2]>60)
dim(newdata)
input<-data[match(rownames(newdata),rownames(data)),]
colMeans(input,na.rm = T)


cor2bed(rownames(input))

rlt<-data.frame(cor2bed(rownames(input)),(rownames(input)))
write.table(rlt,file="4521NASH.candidate.hg19.bed",sep="\t",col.names = F,row.names = F,quote=F)
liftOver 4521NASH.candidate.hg19.bed ~/hpc/db/hg19/hg19ToHg38.over.chain.gz 4521NASH.candidate.hg38.bed unmap




args = commandArgs(trailingOnly=TRUE)
sam<-read.table("~/hpc/epimarker/bedgraph/SampleMatrix.txt",head=T,sep="\t")
data=read.table("pbmc.hg38.matrix.txt",head=T,sep="\t",row.names=1,check.names=F)
filename=unlist(lapply(strsplit(colnames(data),split="[.]"),function(x) x[1]))
newsam=sam[match(filename,sam[,1]),]
blood=which(newsam[,3]=="blood")
solid=which(newsam[,3]=="tissue")
which(rowMeans(data[blood,],na.rm=T)<0.2)




