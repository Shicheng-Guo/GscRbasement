GSIMethDecon<-function(data){
  data<-data.matrix(data)
  group=names(table(colnames(data)))
  index=colnames(data)
  gsi<-gmaxgroup<-avebase<-c()
  for(i in 1:nrow(data)){
  gsit<-0
  gmax<-names(which.max(tapply(as.numeric(data[i,]),index,function(x) mean(x,na.rm=T))))
  for(j in 1:length(group)){
    tmp<-(1-10^(mean(data[i,][which(index==group[j])],na.rm=T))/10^(mean(data[i,][which(index==gmax)],,na.rm=T)))/(length(group)-1)
    gsit<-gsit+tmp
  }
  ave<-tapply(data[i,], index, function(x) mean(x,na.rm=T))
  gmaxgroup<-c(gmaxgroup,gmax)
  gsi<-c(gsi,gsit)
  avebase<-rbind(avebase,ave)
  }
  rlt=data.frame(region=rownames(data),group=gmaxgroup,GSI=gsi,AVE=avebase)
  return(rlt)
}

TopGSIByCategory<-function(gsi,top=1,thresHigh=0.2,thresLow=0.1,allmin=0.05,plotfigure=T,figprefix="tissuespecific"){
  # gsi=GSITest; thresHigh=0.2; thresLow=0.1; allmin=0.05; plotfigure=T; figprefix="test-tissuespecific"
    gsi=data.frame(gsi)
    GSIRlt<-c()
    group<-as.character(unique(gsi$group))
    rank<-c(rep(top,length(group)))
    otf1<-paste(figprefix,"boxplot.pdf",sep="")
    otf2<-paste(figprefix,"heatmap.pdf",sep="")
    parnum<-ceiling(sqrt(length(group)))
    pdf(otf1)
    par(mfrow=c(parnum,parnum),oma = c(2,2,2,2) + 0.1,mar = c(2,2,2,2) + 0.1,cex.axis=0.75, cex.lab=0.75)
    for (i in 1:length(group)){
    # select tissue-specific (remove target group<0.2 or non-target group>0.1)
    subset=gsi[which(gsi$group==group[i]),]
    rexclhigh<-which(apply(subset,1,function(x) x[grep(group[i],colnames(gsi))]<thresHigh))
    rexclall<-which(apply(subset,1,function(x) all(x[grep(group[i],colnames(gsi))]<allmin)))
    xx<-subset[,-grep(group[i],colnames(gsi))]
    rexcllow<-which(apply(xx,1,function(x) any(as.numeric(x[4:length(x)])>thresLow)))
    rexcl<-c(rexclhigh,rexcllow,rexclall)
    subset=subset[-rexcl,]
    subset=subset[order(subset[,3],decreasing=T)[1:rank[i]],]
    GSIRlt<-rbind(GSIRlt,subset)
    if(plotfigure==T){
      zz=subset[which(subset$group==group[i]),]
      if(nrow(zz)>=1){
      boxplot(na.omit(zz[,4:ncol(zz)]),horizontal=T,las=2,col="red")
      }else{
      print(paste(group[i],"do not have biomarker,please check raw data",sep=" "))
      }
    }
  }
  dev.off()
  
  if(plotfigure==T){
    HeatMap(data=data.matrix(na.omit(GSIRlt[,4:ncol(GSIRlt)])),phen=gsub("AVE.","",colnames(GSIRlt)[4:ncol(GSIRlt)]),figure=otf2)
  }
  GSIRlt<-na.omit(GSIRlt)
  return(GSIRlt)
}

