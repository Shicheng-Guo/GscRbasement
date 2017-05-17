
gsi<-function(data){
  group=names(table(colnames(data)))
  index=colnames(data)
  GSI<-c()
  gmaxgroup<-c()
  for(i in 1:nrow(data)){
    gsit<-0
    gmax<-names(which.max(tapply(as.numeric(data[i,]),index,function(x) mean(x,na.rm=T))))
    if(length(gmax)<1){print(data[i,])}
    for(j in 1:length(group)){
      tmp<-(1-10^(mean(na.omit(as.numeric(data[i,which(index==group[j])])),na.rm=T))/10^(mean(na.omit(as.numeric(data[i,which(index==gmax)])))))/(length(group)-1)
      gsit<-c(gsit,tmp)
    }
    gmaxgroup<-c(gmaxgroup,gmax)
    GSI<-c(GSI,sum(gsit,na.rm=T))  # debug for GSI=NAN
  }
  rlt=data.frame(region=rownames(data),group=gmaxgroup,GSI=GSI)
  return(rlt)
}
