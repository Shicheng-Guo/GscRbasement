# EXAMPLE USAGE
#install.packages("gplots")
#install.packages("devtools")
library("gplots")
library("devtools")

sam1<-read.table("GTEx_v7_Annotations_SampleAttributesDS.txt",sep="\t",head=T,fill = T)
sam2<-read.table("GTEx_v7_Annotations_SubjectPhenotypesDS.txt",head=T,sep="\t")
data<-read.table("GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct",sep="\t",head=T,skip=2,row.names=1,check.names = F)
save(data,file="GTEx.RData")
genesymbol<-data[,1]
Data<-data[,-1]
ID=unlist(lapply(strsplit(colnames(Data),"[-]"),function(x) paste(x[1],x[2],sep="-")))


newdata<-Data[match(genesymbol[grep("FGF\\d+$",genesymbol,perl=T)],genesymbol),]
rownames(newdata)<-data[match(rownames(newdata),rownames(data)),1]
sex<-sam2[match(ID,sam2$SUBJID),]$SEX
tissue<-sam1[match(colnames(Data),sam1$SAMPID),]$SMTS
# sampling 7 samples from each catorgy
newdata2<-newdata[,unlist(lapply(unique(as.character(tissue)),function(x) sample(which(tissue %in% x),7)))]
ID=unlist(lapply(strsplit(colnames(newdata2),"[-]"),function(x) paste(x[1],x[2],sep="-")))
sex<-sam2[match(ID,sam2$SUBJID),]$SEX
tissue<-sam1[match(colnames(newdata2),sam1$SAMPID),]$SMTS
newdata3<-t(log(newdata2+1,2))
heatmap3(newdata3)
dev.off()


newdata<-Data[match(genesymbol[grep("FGFR\\d+$",genesymbol,perl=T)],genesymbol),]
rownames(newdata)<-data[match(rownames(newdata),rownames(data)),1]
sex<-sam2[match(ID,sam2$SUBJID),]$SEX
tissue<-sam1[match(colnames(Data),sam1$SAMPID),]$SMTS
# sampling 7 samples from each catorgy
newdata2<-newdata[,unlist(lapply(unique(as.character(tissue)),function(x) sample(which(tissue %in% x),7)))]
ID=unlist(lapply(strsplit(colnames(newdata2),"[-]"),function(x) paste(x[1],x[2],sep="-")))
sex<-sam2[match(ID,sam2$SUBJID),]$SEX
tissue<-sam1[match(colnames(newdata2),sam1$SAMPID),]$SMTS
newdata3<-t(log(newdata2+1,2))
heatmap3(newdata3)
dev.off()


tmp<-read.table("tmp")
newdata<-Data[na.omit(match(tmp[,1],genesymbol)),]
rownames(newdata)<-data[match(rownames(newdata),rownames(data)),1]
sex<-sam2[match(ID,sam2$SUBJID),]$SEX
tissue<-sam1[match(colnames(Data),sam1$SAMPID),]$SMTS
# sampling 7 samples from each catorgy
newdata2<-newdata[,unlist(lapply(unique(as.character(tissue)),function(x) sample(which(tissue %in% x),7)))]
ID=unlist(lapply(strsplit(colnames(newdata2),"[-]"),function(x) paste(x[1],x[2],sep="-")))
sex<-sam2[match(ID,sam2$SUBJID),]$SEX
tissue<-sam1[match(colnames(newdata2),sam1$SAMPID),]$SMTS
newdata3<-t(log(newdata2+1,2))
heatmap3(newdata3)
dev.off()




column_annotation<-data.matrix(data.frame(sex,tissue))
col_annotation<-matrix(data.matrix(data.frame(sex,as.character(tissue))),ncol=2)
# only for FGF1
gene_mRNAseq<-Data[which(data[,1] %in% "FGFR1"),]
data_plot<-tapply(data.matrix(gene_mRNAseq),tissue,function(x) mean(x,na.rm=T))
barplot(data_plot,horiz=T,las=2)
dev.off()

# full panel
row_annotation <- matrix(data.frame(sex,as.character(tissue)))

mat <- matrix(1:50, byrow=T, nrow=5)                   #  x is a matrix of observations by variables.
column_annotation <- matrix(sample(c("red", "blue", "green"), 20, replace=T), ncol=2)
colnames(column_annotation) <- c("Variable X1", "Variable X2")
row_annotation <- matrix(sample(c("red", "blue", "green"), 20, replace=T), nrow=2)
rownames(row_annotation) <- c("Variable Y1", "Variable Y2")
heatmap3(mat, RowSideColors=row_annotation, ColSideColors=column_annotation)
heatmap3(mat,ColSideColors=column_annotation)
