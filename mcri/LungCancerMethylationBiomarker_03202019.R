# DNA methylation research in human lung cancer

## lung cancer (WGBS)
 # 
## lung cancer (RRBS) 
 # 
## lung cancer (450K)
 # GSE39279: 444 lung cancer 450K array
 # GSE92843: 3 lung cancer cell lines (A549, A427 and H322), normal bronchial ephitelial cells (NHBEC)
 # GSE63704: 17 LC,  37 idiopathic lung fibrosis,32 patients suffering from chronic obstructive pulmonary disease and 43 DNA samples derived from healthy-lungs 
 # GSE63940: 36 lung adenocarcinoma cell lines was bisulfite treated and analyzed on the Illumina Infinium HumanMethylation27K 
 # GSE85566: isolated airway epithelial cells of asthmatics and non-asthmatics
## lung cancer (27K)
 # GSE62950: 28 pair LC-Normal(27K and mRNA)
 # GSE32867: 59 pair LC-Normal(27K and mRNA)
 # wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE39nnn/GSE39279/suppl/GSE39279_RAW.tar
 # mkdir GSE39279
 # mv GSE39279_RAW.tar GSE39279
### Illumina GoldenGate Methylation Cancer Panel I
 # GSE16559: 70 normal, 57 lung adenocarcinoma and 158 pleural mesothelioma
 #

library("GEOquery")
# NCI60 with 450K microarray

# 450K with 444 lung cancer 450K array
GSE39279 <- getGEO("GSE39279")
data1 <- as.data.frame(exprs(GSE39279[[1]]))
phen1 <- pData(phenoData(GSE39279[[1]]))
dim(data1)
if(nrow(data1) !=485577) print("data is not complete, check the raw data or idat")

# GSE92843: Three lung cancer cell lines (A549, A427 and H322), normal bronchial ephitelial cells (NHBEC)
GSE92843 <- getGEO("GSE92843")
data2 <- as.data.frame(exprs(GSE92843[[1]]))
phen2 <- pData(phenoData(GSE92843[[1]]))
dim(data2)
if(nrow(data2) !=485577) print("data is not complete, check the raw data or idat")

# GSE63704: 17 LC,  37 idiopathic lung fibrosis,32 patients suffering from chronic obstructive pulmonary disease and 43 DNA samples derived from healthy-lungs 
GSE63704 <- getGEO("GSE63704")
data3 <- as.data.frame(exprs(GSE63704[[1]]))
phen3 <- pData(phenoData(GSE63704[[1]]))
dim(data3)
if(nrow(data3) !=485577) print("data is not complete, check the raw data or idat")

# Illumina GoldenGate Methylation Cancer Panel I
GSE16559 <- getGEO("GSE16559")
data4 <- as.data.frame(exprs(GSE16559[[1]]))
phen4 <- pData(phenoData(GSE16559[[1]]))
dim(data4)
if(nrow(data4) !=485577) print("data is not complete, check the raw data or idat")

# 450K isolated airway epithelial cells of asthmatics and non-asthmatics
GSE85566 <- getGEO("GSE85566")
data5 <- as.data.frame(exprs(GSE85566[[1]]))
phen5 <- pData(phenoData(GSE85566[[1]]))
dim(data5)
if(nrow(data5) !=485577) print("data is not complete, check the raw data or idat")




phen1<-sapply(strsplit(as.character(phen$characteristics_ch1.7),"[:]"),function(x) as.numeric(unlist(x)[2]))  
phen1[phen1==1]<-"Normal"
phen1[phen1==2]<-"schizophrenia"
phen2<-sapply(strsplit(as.character(phen$characteristics_ch1),"[:]"),function(x) (unlist(x)[2]))  # gender
data1=na.omit(data)
PCAPlot(t(data1),phen1,output="GSE41169.scz.normal.pdf",multifigure=T)  # status
PCAPlot(t(data1),phen2,output="GSE41169.gender.pdf",multifigure=T)  # gender



BiocManager::install("MutationalPatterns", version = "3.8")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19", version = "3.8")
library("MutationalPatterns")
library("BSgenome")
head(available.genomes())
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
library(MutationalPatterns)














