# For PCA Analysis to methylation 450K dataset
# for ips methylatin 450K analysis


PCAPlot<-function(data,pheno,output,multifigure=T){
  pca <- prcomp(data,center=T,scale = F)  # Here, input file: row is individual and column is variable
  outputfile=paste(output,".pdf",sep="")
  pdf(outputfile)
  if(multifigure){
    par(mfrow=c(2,2),mar=c(4,4,4,4))  
  }
  plot((pca$sdev[1:10])^2,type="o",xaxt="n",ylab="Variances",xlab="Principle Components",col="red",lwd=2)
  axis(1,at=0:10,labels=paste("PC",0:10,sep=""))
  var<-c()
  for(i in 1:length(pca$sdev)){var[i]<-sum((pca$sdev[1:i])^2)/sum((pca$sdev)^2)}
  plot(var,ylab="total variance",xlab="number of principle components",lwd=2,type="l")
  abline(h=0.8,col="grey",lty=2)
  abline(v=which(var>0.8)[1],col="grey",lty=2)
  scores <- data.frame(pheno, pca$x[,1:3])
  col = as.numeric(as.factor(pheno))
  plot(x=scores$PC1,y=scores$PC2, xlim=c(min(scores$PC1),max(scores$PC1)),ylim=c(min(scores$PC2),max(scores$PC2)),type="n",xlab="PC1",ylab="PC2")
  for(i in 1:length(scores$PC1)){
    points(scores$PC1[i],scores$PC2[i],pch=as.numeric(as.factor(pheno))[i],col=col[i],cex=0.8,lwd=2)
  }
  legend("bottomright",legend=names(table(pheno)),pch=1:length(table(pheno)),col=1:length(table(pheno)),bty="n")
  plot(x=scores$PC1,y=scores$PC3, xlim=c(min(scores$PC1),max(scores$PC1)),ylim=c(min(scores$PC3),max(scores$PC3)),type="n",xlab="PC1",ylab="PC3")
  for(i in 1:length(scores$PC1)){
    points(scores$PC1[i],scores$PC3[i],pch=as.numeric(as.factor(pheno))[i],col=col[i],cex=0.9,lwd=2)
  }
  legend("bottomright",legend=names(table(pheno)),pch=1:length(table(pheno)),col=1:length(table(pheno)),bty="n")
  dev.off()
}


library("GEOquery")
Access="GSE82218"
data <- getGEO(Access,destdir="/home/shg047/NAS2/GEO")
Rawsave<-paste(Access,"_matrix.Rdata",sep="")
save(data, file=Rawsave)

load(Rawsave)
beta <- as.data.frame(exprs(data[[1]]))
phen <- pData(phenoData(data[[1]]))

phen1<-unlist(lapply(as.character(phen$characteristics_ch1),function(x) unlist(strsplit(x,"[: ]"))[length(unlist(strsplit(x,"[: ]")))]))
phen2<-unlist(lapply(as.character(phen$characteristics_ch1.1),function(x) unlist(strsplit(x,"[: ]"))[length(unlist(strsplit(x,"[: ]")))]))
phen3<-unlist(lapply(as.character(phen$characteristics_ch1.2),function(x) unlist(strsplit(x,"[: ]"))[length(unlist(strsplit(x,"[: ]")))]))

# assess the imputation accuracy
library("impute")
BETA<-as.numeric(data.matrix(beta))
Sample<-sample(1:length(BETA),2000)
Value<-BETA[Sample]
BETA[Sample]<-NA
BETA2<-matrix(BETA,ncol=ncol(beta),byrow=F)
BETA.normal<-impute.knn(data.matrix(BETA2),rowmax = 0.8)
BETA.normal2<-as.numeric(BETA.normal$data)
Value2<-BETA.normal2[Sample]
Error<-Value-Value2
pdf("impute.error.pdf")
plot(hist(Error),col="blue")
dev.off()

beta.normal<-impute.knn(data.matrix(beta),rowmax = 0.8)
Beta=na.omit(beta.normal$data)
PCAPlot(t(Beta),phen1,output="PCA.phen1.pdf",multifigure=T)
PCAPlot(t(Beta),phen2,output="PCA.phen2.pdf",multifigure=T)
PCAPlot(t(Beta),phen3,output="PCA.phen3.pdf",multifigure=T)

