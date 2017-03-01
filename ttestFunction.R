ttestFunction<-function(data,x1,x2){
  output<-matrix(NA,dim(data)[1],4)
  data=data+matrix(abs(rnorm(nrow(data)*ncol(data),0,0.0001)),nrow(data),ncol(data))
  for(i in 1:dim(data)[1]){
    out<-data.frame()
    if(all(! any(all(is.na(data[i,x1])),sum(! is.na(data[i,x1]))<2,sum(! is.na(data[i,x2]))<2,all(is.na(data[i,x2]))),sum(is.na(data[i,]))<0.5*length(data[i,]))){ 
      tmp1<-try(t.test(as.numeric(data[i,x1]),as.numeric(data[i,x2]), na.action=na.omit))
      output[i,1]<-tmp1$p.value
      output[i,2]<-as.numeric((mean(as.numeric(data[i,x1]),na.rm=T)-mean(as.numeric(data[i,x2]),na.rm=T)))
      output[i,3]<-mean(as.numeric(data[i,x1]),na.rm=T)
      output[i,4]<-mean(as.numeric(data[i,x2]),na.rm=T)
      print(i)
    }
  }
  rownames(output)<-rownames(data)
  P.Adj<-p.adjust(output[,1],method="fdr")
  out<-data.frame(output[,1],P.Adj,output[,2:4])
  colnames(out)=c("P-value","P.FDR","Delta","G1","G2")
  return(out)
}
