setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation")
library("easyPubMed")
Count<-c()
Pubid<-c()
RAG<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/hg19/refGene.hg19.V2.bed",head=F,sep="\t")
RAG<-unique(RAG[,5])

for(i in 1:length(RAG)){
  gene<-RAG[i]
  my_query <- paste("methylation[TIAB] AND meta[TIAB] AND",gene,"[TIAB]")
  my_query <- get_pubmed_ids(pubmed_query_string = my_query)
  Count<-rbind(Count,c(as.character(RAG[i]),my_query$Count,paste(unlist(my_query$IdList),collapse=";")))
  print(c(as.character(RAG[i]),my_query$Count))
}
write.table(Count,file="methylation.meta.txt",sep="\t",quote=F)
