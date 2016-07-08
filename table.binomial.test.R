table.binomial.test<-function(region,event.time){
  region1<-table(region)
  region2<-unique(region)
  p=length(region)/(length(region2)*event.time)
  pvalue<-c()
  for(i in 1:length(table(region))){
    ptmp<-binom.test(table(region)[i],event.time,p,alternative="two.sided")$p.value
    pvalue<-c(pvalue,ptmp)
  }
  return(pvalue)
  print(paste("there are",sum(pvalue<0.05/length(pvalue),"regions were over-preferred in the samples"),sep=" "))
}
