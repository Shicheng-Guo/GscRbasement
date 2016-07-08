cor2map<-function(cor,extend=1){
  a<-unlist(lapply(strsplit(as.character(cor),split=c(":")),function(x) x))
  bed<-matrix(a,ncol=2,byrow=T)
  bed<-data.frame(start=bed,end=as.numeric(bed[,2])+extend)
  return(bed)
}
