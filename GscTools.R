
ENST2Symbol<-function(ENST){
  db<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/ENSG.ENST.ENSP.Symbol.hg19.bed",sep="\t")
  Symbol<-db[match(ENST,db$V7),4]
  return(Symbol)
}

ENSG2Symbol<-function(ENSG){
  db<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/ENSG.ENST.ENSP.Symbol.hg19.bed",sep="\t")
  Symbol<-db[match(ENSG,db$V8),4]
  return(Symbol)
}

ENSP2Symbol<-function(ENSP){
  db<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/ENSG.ENST.ENSP.Symbol.hg19.bed",sep="\t")
  Symbol<-db[match(ENSP,db$V9),4]
  return(Symbol)
}
