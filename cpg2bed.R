cpg2bed<-function(cpg){
map<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/hg19/GPL13534_450K_hg19_V3.bed")
bed<-map[match(cpg,map[,4]),]
return(bed)
}
