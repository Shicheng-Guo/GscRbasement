PairTtestPValue<-function(data,x1,x2,pair=FALSE){
  data<-data.matrix(data)
  output<-matrix(NA,dim(data)[1],6)   # set output matrix ()
  for(i in 1:dim(data)[1]){
    out<-data.frame()
    if(pair==TRUE){
      Valid<-nrow(na.omit(data.frame(data[i,x1],data[i,x2])))
    }else{
      Valid<-100
    }
    if( sum(!is.na(data[i,x1]))>=3 & sum(!is.na(data[i,x2]))>=3 & Valid>3){ 
      tmp1<-try(t.test((data[i,x1]),(data[i,x2]),paired=pair, na.action=na.omit))
      output[i,1]<-format(tmp1$p.value, scientific=TRUE)
      output[i,2]<-round(mean((data[i,x1]))-mean((data[i,x2])),3)
      output[i,3]<-round(mean((data[i,x1])),3)
      output[i,4]<-round(mean((data[i,x2])),3)
      output[i,5]<-round(sd(data[i,x1]),3)
      output[i,6]<-round(sd(data[i,x2]),3)
      print(i)
    }
  }
  rownames(output)<-rownames(data)
  output
}
