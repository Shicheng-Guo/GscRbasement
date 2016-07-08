corsort<-function(cor){
  a<-unlist(lapply(strsplit(as.character(cor),split=c(":")),function(x) x))
  bed<-matrix(a,ncol=2,byrow=T)
  bed<-bed[order(bed[,1],as.numeric(bed[,2])),]
  cor<-apply(bed,1,function(x){paste(unlist(strsplit(x,"\t"))[1],":",unlist(strsplit(x,"\t"))[2],sep="")})
  return(cor)
}
