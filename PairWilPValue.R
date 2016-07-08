PairWilPValue<-function(data,x1,x2){
  output<-matrix(NA,dim(data)[1],5)
  for(i in 1:dim(data)[1]){
    out<-data.frame()
    if(all(! any(all(is.na(data[i,x1])),all(is.na(data[i,x2]))),sum(is.na(data[i,]))<0.5*length(data[i,]))){ 
      tmp1<-try(wilcox.test(as.numeric(data[i,x1]),as.numeric(data[i,x2]),paired=T, na.action=na.omit))
      output[i,1]<-tmp1$p.value
      output [i,2]<-mean(as.numeric(data[i,x1]))-mean(as.numeric(data[i,x2]))
      output[i,3]<-"wilcox"
      output[i,4]<-mean(as.numeric(data[i,x1]))
      output[i,5]<-mean(as.numeric(data[i,x2]))
      # print(i)
    }
  }
  out<-cbind(rownames(data),output)
  out
}
