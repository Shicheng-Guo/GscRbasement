library(Biobase)
library(GEOquery)
library(limma)

GSE125362 <- getGEO("GSE125362", GSEMatrix =TRUE, AnnotGPL=TRUE)
GSE125362 <- as.data.frame(exprs(GSE125362[[1]]))
phen <- pData(phenoData(GSE125362[[1]]))
head(GSE125362[1:5,1:5])

GSE32413 <- getGEO("GSE32413", GSEMatrix =TRUE, AnnotGPL=TRUE)
GSE32413 <- as.data.frame(exprs(GSE32413[[1]]))
phen <- pData(phenoData(GSE32413[[1]]))
head(GSE32413[1:5,1:5])

GSE125362 <- getGEO("GSE125362", GSEMatrix =TRUE, AnnotGPL=TRUE)
GSE125362 <- as.data.frame(exprs(GSE125362[[1]]))
phen <- pData(phenoData(GSE125362[[1]]))
head(GSE125362[1:5,1:5])

GSE76885 <- getGEO("GSE76885", GSEMatrix =TRUE, AnnotGPL=TRUE)
GSE76885 <- as.data.frame(exprs(GSE76885[[1]]))
phen <- pData(phenoData(GSE76885[[1]]))
head(GSE76885[1:5,1:5])

GSE76807 <- getGEO("GSE76807", GSEMatrix =TRUE, AnnotGPL=TRUE)
GSE76807 <- as.data.frame(exprs(GSE76807[[1]]))
phen <- pData(phenoData(GSE76807[[1]]))
head(GSE76807[1:5,1:5])

GSE68698 <- getGEO("GSE68698", GSEMatrix =TRUE, AnnotGPL=TRUE)
GSE68698 <- as.data.frame(exprs(GSE68698[[1]]))
phen <- pData(phenoData(GSE68698[[1]]))
head(GSE68698[1:5,1:5])

GSE59785 <- getGEO("GSE59785", GSEMatrix =TRUE, AnnotGPL=TRUE)
GSE59785 <- as.data.frame(exprs(GSE59785[[1]]))
phen <- pData(phenoData(GSE59785[[1]]))
head(GSE59785[1:5,1:5])

GSE45485 <- getGEO("GSE45485", GSEMatrix =TRUE, AnnotGPL=TRUE)
GSE45485 <- as.data.frame(exprs(GSE45485[[1]]))
phen <- pData(phenoData(GSE45485[[1]]))
head(GSE45485[1:5,1:5])

GSE19617 <- getGEO("GSE19617", GSEMatrix =TRUE, AnnotGPL=TRUE)
GSE19617 <- as.data.frame(exprs(GSE19617[[1]]))
phen <- pData(phenoData(GSE19617[[1]]))
head(GSE19617[1:5,1:5])

GSE76809 <- getGEO("GSE19617", GSEMatrix =TRUE, AnnotGPL=TRUE)
GSE76809 <- as.data.frame(exprs(GSE76809[[1]]))
phen <- pData(phenoData(GSE76809[[1]]))
head(GSE76809[1:5,1:5])



dim(GSE125362)
dim(GSE32413)
dim(GSE76885)
dim(GSE76807)
dim(GSE68698)
dim(GSE59785)
dim(GSE45485)
dim(GSE19617)

system("wget https://cran.r-project.org/src/contrib/ff_2.2-14.tar.gz")
system("wget https://cran.r-project.org/src/contrib/ffbase_0.12.7.tar.gz")
system("wget https://cran.r-project.org/src/contrib/SuppDists_1.1-9.4.tar.gz")

install.packages("ff_2.2-14.tar.gz")
install.packages("ffbase_0.12.7.tar.gz")
install.packages("SuppDists_1.1-9.4.tar.gz")
fastmatch


library("Haplin")

data1<-read.table("Iron.Gene.in.SSc.Trim.Pvalue.txt",head=T,sep="\t")
# house-keeping
hkg<-read.table("Human_Housekeeping_Gene_List.txt")
match(hkg[,1],)
png("qqplot.png")
pQQ(data$P.Value, nlabs =nrow(data), conf = 0.95)
dev.off()

pQQ


phen1<-sapply(strsplit(as.character(phen$characteristics_ch1.7),"[:]"),function(x) as.numeric(unlist(x)[2]))  # status 1:control, 2:scz
phen1[phen1==1]<-"Normal"
phen1[phen1==2]<-"schizophrenia"
phen2<-sapply(strsplit(as.character(phen$characteristics_ch1),"[:]"),function(x) (unlist(x)[2]))  # gender

data1=na.omit(data)
PCAPlot(t(data1),phen1,output="GSE41169.scz.normal.pdf",multifigure=T)  # status
PCAPlot(t(data1),phen2,output="GSE41169.gender.pdf",multifigure=T)  # gender


