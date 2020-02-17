BiocManager::install("KEGGprofile")
library(limma)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(KEGGprofile)

data(pro_pho_expr)
data(pho_sites_count)
 
tab <- getGeneKEGGLinks(species="hsa")
xtab<- mapIds(org.Hs.eg.db, tab$GeneID,column="SYMBOL", keytype="ENTREZID")
head(tab)

setwd("/mnt/bigdata/Genetic/Projects/shg047/meta/panc754")
download_KEGGfile(pathway_id="04110",species='hsa')
download_KEGGfile(pathway_id="05323",species='hsa')

genes<-row.names(pho_sites_count)[which(pho_sites_count>=10)]
pho_KEGGresult<-find_enriched_pathway(genes,species='hsa')
pho_KEGGresult[[1]][,c(1,5)]

pho_expr<-pro_pho_expr[,7:12]
temp<-apply(pho_expr,1,function(x) length(which(is.na(x))))
pho_expr<-pho_expr[which(temp==0),]
col<-col_by_value(pho_expr,col=colorRampPalette(c('green','black','red'))(1024),range=c(-6,6))
temp<-plot_pathway(pho_expr,type="bg",bg_col=col,text_col="white",magnify=1.2,species='hsa',database_dir=system.file("extdata",package="KEGGprofile"),pathway_id="04110")

col<-col_by_value(pho_sites_count,col=colorRampPalette(c('white','khaki2'))(4),breaks=c(0,1,4,10,Inf)) ## visualization by method 'lines'
temp<-plot_pathway(pro_pho_expr,type="lines",bg_col=col,line_col=c("brown1","seagreen3"),groups=c(rep("Proteome",6),rep("Phosphoproteome",6)),magnify=1.2,species='hsa',database_dir=system.file("extdata",package="KEGGprofile"),pathway_id="04110",max_dist=5)

