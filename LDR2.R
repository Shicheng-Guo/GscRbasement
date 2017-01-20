
library("genetics")
hapinfo<-c(rep("CCCC",8),rep("TTTT",2))

# average linkage among all CpG locis in haplotype
hapinfo2LDR2(hapinfo)

# circle plot base on hapinfo 
bspplot(hapinfo2matrix(hapinfo))

# linkage of loci 1 and loci 2
diseq(genotype(metype2Geno(hapinfo)[[1]]))

metype2Geno<-function(hapinfo){
rlt<-list()
for(i in 1:(nchar(hapinfo)[1]-1)){
  Metype<-unlist(lapply(hapinfo,function(x) substr(x,i,i+1)))
  Geno<-lapply(Metype,meth2Geno)
  rlt[[i]]<-as.vector(unlist(Geno))
}
return(rlt)
}

hapinfo2LDR2<-function(hapinfo){
  # R2
  # Dpprime to indicate strength of CC and TT, not same as genetic define (CT)
  rlt<-list()
  Mhap<-lapply(metype2Geno(hapinfo),function(x){diseq(genotype(as.vector(x)))})
  R2<-mean(unlist(lapply(Mhap,function(x) x$R2.overall)),na.rm=T)
  Dpprime<-mean(unlist(lapply(Mhap,function(x) x$Dprime.overall)),na.rm=T)
  rlt$R2<-R2
  rlt$Dpprime<--Dpprime  
  return(rlt)
}

meth2Geno<-function(x){
  as.list(paste(as.character(unlist(strsplit(x,split=""))),collapse="/"))
}

hapinfo2matrix<-function(x){
  matrix<-do.call("rbind",as.list(hapinfo))
  matrix<-data.matrix(do.call("rbind",sapply(matrix,function(x) strsplit(x,""))))
  matrix[matrix=="C"]=1
  matrix[matrix=="T"]=0
  rownames(matrix)<-make.names(rownames(matrix),unique = TRUE)
  matrix<-matrix(as.numeric(matrix),nrow=nrow(matrix),byrow=F)
  matrix
}

bspplot<-function(Matrix){
  par(mar=c(3,3,1,5))
  col=colorRampPalette(c("white", "red"))(20)
  circle=c(1,19)
  plot(x=nrow(Matrix),y=ncol(Matrix),type="n",xlab="",ylab="",xlim=c(0,ncol(Matrix)+1),ylim=c(0,nrow(Matrix)+1))
  for(i in 1:ncol(Matrix)){
    for(j in 1:nrow(Matrix)){
      points(i,j,col=1,pch=circle[Matrix[j,i]+1],cex=1)
    }
  }
}



