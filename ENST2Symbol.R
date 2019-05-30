ENST2Symbol<-function(ENST){
  db<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/ENSG.ENST.ENSP.Symbol.hg19.bed",sep="\t")
  Symbol<-db[match(ENST,db$V7),4]
  return(Symbol)
}
