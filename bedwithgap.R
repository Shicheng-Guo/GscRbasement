bedwithgap<-function(bed,gap){
  bed<-as.matrix(bed)
  bed[,2]=as.numeric(bed[,2])-gap
  bed[,3]=as.numeric(bed[,3])+gap
  bed<-data.frame(bed)
  bed
}
