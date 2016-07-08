ttestFunction<-function(data,x1,x2){
  output<-matrix(NA,dim(data)[1],5)
  for(i in 1:dim(data)[1]){
    out<-data.frame()
    if(all(! any(all(is.na(data[i,x1])),length(na.omit(data[i,x1]))<2,length(na.omit(data[i,x2]))<2,all(is.na(data[i,x2]))),sum(is.na(data[i,]))<0.5*length(data[i,]))){ 
      tmp1<-try(t.test(as.numeric(data[i,x1]),as.numeric(data[i,x2]), na.action=na.omit))
      output[i,1]<-tmp1$p.value
      output[i,2]<-as.numeric((mean(as.numeric(data[i,x1]),na.rm=T)-mean(as.numeric(data[i,x2]),na.rm=T)))
      output[i,3]<-"t-test"
      output[i,4]<-mean(as.numeric(data[i,x1]),na.rm=T)
      output[i,5]<-mean(as.numeric(data[i,x2]),na.rm=T)
      print(i)
    }
  }
  out<-cbind(rownames(data),output)
  P.Adj<-p.adjust(output[,1],method="fdr")
  out<-data.frame(out[,1],out[,2],P.Adj,out[,3:6])
  colnames(out)=c("RowName","P-value","P.FDR","Delta","Test","M(G1)","M(G2)")
  return(out)
}
