AverageWithGroup<-function(data){
  base=unlist(lapply(strsplit(colnames(data),"[.]"),function(x) x[[1]]))
  matrix=apply(data,1,function(x) tapply(x,base,function(x) mean(x,na.rm=T)))
  matrix<-t(matrix)
  rownames(matrix)=rownames(data)
  matrix<-matrix[!rowSums(!is.finite(matrix)),]
  return(matrix)
}
