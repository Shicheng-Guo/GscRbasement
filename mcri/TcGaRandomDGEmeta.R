
setwd("/mnt/bigdata/Genetic/Projects/shg047/temp/housekeeping")
data<-read.csv("https://raw.githubusercontent.com/Shicheng-Guo/esophageal/master/extdata/ESCC-UKBB-RSKAT2020.csv",as.is=T)
data<-read.csv("https://raw.githubusercontent.com/Shicheng-Guo/colon/master/extdata/UKBB50K/CRC-UKBB-RSKAT2020.csv",as.is=T)
data<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/housekeeping/master/housekeeping.txt",as.is=T,head=T,sep="\t")
symbol<-data[,1]
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/PancancerMetaDge.R")
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/pancancermetaOsHr.R")
pancancermetadge(symbol,"housekeeping")
pancancermetaOsHr(symbol,"housekeeping")
OROSmerge("housekeep")

setwd("~/hpc/methylation/chol/panmeta")
memo="cholangiocarcinoma"
data<-read.csv("https://raw.githubusercontent.com/Shicheng-Guo/esophageal/master/extdata/ESCC-UKBB-RSKAT2020.csv",as.is=T)
data<-read.csv("https://raw.githubusercontent.com/Shicheng-Guo/colon/master/extdata/UKBB50K/CRC-UKBB-RSKAT2020.csv",as.is=T)
data<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/housekeeping/master/housekeeping.txt",as.is=T,head=T,sep="\t")
data<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/cholangiocarcinoma/master/cholTargetRegion50.txt",as.is=T)
symbol<-data[,1]
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/PancancerMetaDge.R")
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/pancancermetaOSHR.R")
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/OROSmerge.R")
pancancermetadge(symbol,memo)
pancancermetaOSHR(symbol,memo)
OROSmerge(memo)

setwd("/home/guosa/hpc/methylation/brca/meta")
memo="breastcancer"
data<-read.csv("https://raw.githubusercontent.com/Shicheng-Guo/esophageal/master/extdata/ESCC-UKBB-RSKAT2020.csv",as.is=T)
data<-read.csv("https://raw.githubusercontent.com/Shicheng-Guo/colon/master/extdata/UKBB50K/CRC-UKBB-RSKAT2020.csv",as.is=T)
data<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/housekeeping/master/housekeeping.txt",as.is=T,head=T,sep="\t")
data<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/cholangiocarcinoma/master/cholTargetRegion50.txt",as.is=T)
data<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/brca/master/extdata/brca.tcga.target.hg19.bed",as.is=T)
head(data)
symbol<-data[,4]
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/Pancancermetadge.R")
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/pancancermetaOsHr.R")
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/OROSmerge.R")
pancancermetadge(symbol,memo)
pancancermetaOSHR(symbol,memo)
OROSmerge(memo)







cpg<-data[grep("cg",data[,1]),]
cpg2symbol(as.character(cpg))
bed<-cpg2bed(as.character(cpg))

cpg2symbol<-function(cpg){
map<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/hg19/GPL13534_450K_hg19_V3.bed")
symbol<-map[which(map[,4] %in% cpg),5]
return(symbol)
}

cpg2bed<-function(cpg){
map<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/hg19/GPL13534_450K_hg19_V3.bed")
bed<-map[match(cpg,map[,4]),]
return(bed)
}
