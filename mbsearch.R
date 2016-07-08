mbsearch<-function(data,sortedmapfile,window=100,mincpgnumber=4){
  # rowname of data is cpg site and column name is sample id
  # map file is bed file and the fourth column is cpg 
  output<-list()
  data<-RawNARemove(data)
  data<-impute.knn(data.matrix(data))$data
  map<-sortedmapfile
  map<-map[map[,4] %in% rownames(data),]
  newdata<-data[match(map[,4],rownames(data)),]
  if(nrow(map)==nrow(newdata)){
    a<-map[,2]
    i=1
    tmp<-c()
    rlt<-c()
    rowname<-c()
    index<-0
    while(i < length(a)){
      start=a[i]
      end=a[i+1]
      end-start
      if(end-start<window && end-start>0){
        tmp<-c(tmp,i)
        i=i+1
      }else{
        if(length(tmp)>mincpgnumber){
          index=index+1
          tmp<-c(min(tmp),max(tmp),length(tmp),a[max(tmp)]-a[min(tmp)],round(length(tmp)/(a[max(tmp)]-a[min(tmp)]),4),round((a[max(tmp)]-a[min(tmp)])/(length(tmp))))
          rlt<-rbind(rlt,tmp)
          tmp2<-paste(map[tmp[1],1],":",map[tmp[1],2],"-",map[tmp[2],2],sep="")
          rowname<-c(rowname,tmp2)
        }
        tmp<-c()
        i=i+1
      }
    }
    rownames(rlt)<-rowname
    colnames(rlt)<-c("rowstart","rowstart","CpGNumber","RegionLength","CpGRatio","AverageGap")
    cor<-c()
    for(j in 1:nrow(rlt)){
      cor1<-mean(cor(t(newdata[rlt[j,1]:rlt[j,2],seq(1,ncol(newdata),by=2)]),use="complete.obs")) # cancer
      cor2<-mean(cor(t(newdata[rlt[j,1]:rlt[j,2],seq(2,ncol(newdata),by=2)]),use="complete.obs")) # normal
      tmp<-c(cor1,cor2)
      cor<-rbind(cor,tmp)
    }
    rownames(cor)<-rowname
    colnames(cor)<-c("ratioA","ratioB")
    
    output$hdr<-rlt
    output$hdrc<-cor
    
    return(output)
  }
}
