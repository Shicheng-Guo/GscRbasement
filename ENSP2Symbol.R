ENSP2Symbol<-function(ENSP){
  db<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/ENSG.ENST.ENSP.Symbol.hg19.bed",sep="\t")
  Symbol<-db[match(ENSP,db$V9),4]
  return(Symbol)
}
