map2cor<-function(map){
  cor<-apply(map,1,function(x){paste(unlist(strsplit(x,"\t"))[1],":",unlist(strsplit(x,"\t"))[2],sep="")})
  cor<-gsub("[ ]","",cor)
  return(cor)
}
