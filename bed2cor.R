
bed2cor<-function(bed){
cor<-apply(bed,1,function(x) paste(x[1],":",x[2],"-",x[3],sep=""))
return(cor)
}
