library("easyPubMed")
Count<-c()
Pubid<-c()
RAG<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/hg19/refGene.hg19.V2.bed",head=F,sep="\t")
RAG<-unique(RAG[,5])
for(i in 1:length(RAG)){
gene<-RAG[i]
my_query <- paste("(SNP[TIAB] or polymorphism[TIAB]) AND association[TIAB] AND rheumatoid arthritis[TIAB] AND",gene,"[TIAB]")
my_query <- get_pubmed_ids(pubmed_query_string = my_query)
Count<-rbind(Count,c(as.character(RAG[i]),my_query$Count,paste(unlist(my_query$IdList),collapse=";")))
print(c(as.character(RAG[i]),my_query$Count))
}
write.table(Count,file="RA.polymorphism.meta.txt",sep="\t",quote=F)
