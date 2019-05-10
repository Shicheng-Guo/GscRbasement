# the easilest cross-validation to glm prediction
library("MASS")
library("readxl")
library("boot")
setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/HCC")
source("https://raw.githubusercontent.com/Shicheng-Guo/Gscutility/master/assess.R")
beta<-read.table("methyLevel_54samples.txt",head=T,sep="\t")
phen<-data.frame(read_xlsx("ROC.xlsx",sheet=1))
rownames(beta)<-paste(beta[,1],":",beta[,2],"-",beta[,3],sep="")
beta<-beta[,4:ncol(beta)]
input<-na.omit(data.frame(phen=phen$phen,t(beta)))
set.seed(200)
cost <- function(r, pi = 0) mean(abs(r-pi) > 0.5)
glm.fit <- glm(phen~.,data=input)
xx<-cv.glm(input,glm.fit,cost,K=10)$delta
