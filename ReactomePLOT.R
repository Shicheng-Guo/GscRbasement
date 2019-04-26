ReactomePLOT<-function(symbols){
library("ReactomePA")
library("org.Hs.eg.db")
ENTREZID<-mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')
x <- enrichPathway(gene=ENTREZID,pvalueCutoff=0.05, readable=T)
barplot(x, showCategory=30)
dotplot(x, showCategory=30)
emapplot(x)
cnetplot(x, categorySize="pvalue", foldChange=geneList)
}
