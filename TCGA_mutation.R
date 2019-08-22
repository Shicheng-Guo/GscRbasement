BiocManager::install("TCGAbiolinks")
BiocManager::install("DT")

library("TCGAbiolinks")
library("DT")

maf <- GDCquery_Maf("CHOL", pipelines = "muse")
maf[1:20,]
datatable(maf[1:20,],
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

query.maf.hg19 <- GDCquery(project = "TCGA-CHOL", 
                           data.category = "Simple nucleotide variation", 
                           data.type = "Simple somatic mutation",
                           access = "open", 
                           legacy = TRUE)
library(maftools)
library(dplyr)
query <- GDCquery_Maf("CHOL", pipelines = "muse")

newdata<-data.frame(query$Tumor_Sample_Barcode,query$Hugo_Symbol,query$Variant_Classification,query$Variant_Type)
head(newdata)

maf <- GDCquery_Maf("CHOL", pipelines = "muse") %>% read.maf
datatable(getSampleSummary(maf),
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)
plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
oncoplot(maf = maf, top = 10, removeNonMutated = TRUE)
titv = titv(maf = maf, plot = FALSE, useSyn = TRUE)
plotTiTv(res = titv)








