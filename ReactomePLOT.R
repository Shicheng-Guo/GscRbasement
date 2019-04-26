ReactomePLOT<-function(symbols,output){
library("ReactomePA")
library("org.Hs.eg.db")
ENTREZID<-mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')
x <- enrichPathway(gene=ENTREZID,pvalueCutoff=0.05, readable=T)
P1<-barplot(x, showCategory=30)
P2<-dotplot(x, showCategory=30)
P3<-emapplot(x)
ggsave(P1,file=paste(output,"reactome.barplot.pdf",sep="."))
ggsave(P2,file=paste(output,"reactome.dotplot.pdf",sep="."))
ggsave(P3,file=paste(output,"reactome.emapplot.pdf",sep="."))
}
