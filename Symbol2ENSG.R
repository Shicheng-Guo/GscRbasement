Symbol2ENSG<-function(Symbol){
  db<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/ENSG.ENST.ENSP.Symbol.hg19.bed",sep="\t")
  ENSG<-na.omit(as.character(db[match(Symbol,db$V4),8]))
  return(ENSG)
}

# vip<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/rheumatology/pharmacogenomics/RA/VIP.Gene.hg19.bed",head=F)
# Symbol2ENSG(vip$V4)
