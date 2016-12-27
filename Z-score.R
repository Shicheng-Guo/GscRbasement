setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod\\analysis\\tissue-of-origin-mapping")
normal<-data.matrix(read.table("normal-ref.txt",sep="\t",head=T,row.names=1,as.is=T))
CCP<-data.matrix(read.table("CCP.txt",sep="\t",head=T,row.names=1,as.is=T))
LCP<-data.matrix(read.table("LCP.txt",sep="\t",head=T,row.names=1,as.is=T))

#################################################
################## CCP ##########################
#################################################
Zmax<-matrix(nrow=nrow(CCP),ncol=ncol(CCP))
for(i in 1:ncol(CCP)){
  for(k in 1:nrow(normal)){
    idx<-1:ncol(normal)
    idx1<-sample(idx,30)
    idx2<-idx[which(! idx %in% idx1)]
    
    Mean<-mean(normal[k, idx1])
    SD<-sd(normal[k, idx1])
    
    Z<-c()
    for(p in normal[k,idx2]){
      z <- (p - Mean)/(SD/sqrt(length(idx1)))
      Z<-c(Z,z)
    }
    zmp <- (CCP[k,i] - mean(normal[k, idx]))/(sd(normal[k, idx])*sqrt((length(normal[k, idx])-1)/(length(normal[k, idx]))))
    Zmax[k,i]=zmp 
  }
}
rownames(Zmax)=rownames(CCP)
colnames(Zmax)=colnames(CCP)
Zmax

par(mfrow=c(3,1))
library("beeswarm")
input<-list(Brain=Zmax[1,],Colon=Zmax[2,],Intestine=Zmax[3,],Kidney=Zmax[4,],Liver=Zmax[5,],Lung=Zmax[6,],Pancreas=Zmax[7,],Spleen=Zmax[8,],stomach=Zmax[9,],WBC=Zmax[10,],CT=Zmax[11,])
beeswarm(input,col = 2:12,pch=16,method="center",ylab="Z-Score",main="CRC plamsa")


Zmax<-matrix(nrow=nrow(LCP),ncol=ncol(LCP))
for(i in 1:ncol(LCP)){
  for(k in 1:nrow(normal)){
    idx<-1:ncol(normal)
    idx1<-sample(idx,29)
    idx2<-idx[which(! idx %in% idx1)]
    
    Mean<-mean(normal[k, idx1])
    SD<-sd(normal[k, idx1])
    
    Z<-c()
    for(p in normal[k,idx2]){
      z <- (p - Mean)/(SD/sqrt(length(idx1)))
      Z<-c(Z,z)
    }
    
    zmp <- (LCP[k,i] - mean(normal[k, idx]))/(sd(normal[k, idx])*sqrt((length(normal[k, idx])-1)/(length(normal[k, idx]))))
    Zmax[k,i]=zmp 
  }
}

rownames(Zmax)=rownames(LCP)
colnames(Zmax)=colnames(LCP)
Zmax

library("beeswarm")
input<-list(Brain=Zmax[1,],Colon=Zmax[2,],Intestine=Zmax[3,],Kidney=Zmax[4,],Liver=Zmax[5,],Lung=Zmax[6,],Pancreas=Zmax[7,],Spleen=Zmax[8,],stomach=Zmax[9,],WBC=Zmax[10,],CT=Zmax[11,])
beeswarm(input,col = 2:12,pch=16,method="center",ylab="Z-Score",main="LC plamsa")


Zmax<-matrix(nrow=nrow(normal),ncol=ncol(normal))
for(i in 1:ncol(normal)){
  for(k in 1:nrow(normal)){
    idx<-1:ncol(normal)
    idx1<-sample(idx,75)
    idx2<-idx[which(! idx %in% idx1)]
    
    Mean<-mean(normal[k, idx1])
    SD<-sd(normal[k, idx1])
    
    Z<-c()
    for(p in normal[k,idx2]){
      z <- (p - Mean)/(SD/sqrt(length(idx1)))
      Z<-c(Z,z)
    }
    
    zmp <- (normal[k,i] - mean(normal[k, idx]))/(sd(normal[k, idx])*sqrt((length(normal[k, idx])-1)/(length(normal[k, idx]))))
    Zmax[k,i]=zmp 
  }
}

rownames(Zmax)=rownames(normal)
colnames(Zmax)=colnames(normal)
Zmax

library("beeswarm")
input<-list(Brain=Zmax[1,],Colon=Zmax[2,],Intestine=Zmax[3,],Kidney=Zmax[4,],Liver=Zmax[5,],Lung=Zmax[6,],Pancreas=Zmax[7,],Spleen=Zmax[8,],stomach=Zmax[9,],WBC=Zmax[10,],CT=Zmax[11,])
beeswarm(input,col = 2:12,pch=16,method="center",ylab="Z-Score",main="Normal plamsa")

