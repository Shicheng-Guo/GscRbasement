#############################
#     LAMBDA ESTIMATION     #
#############################

# Using the library “gap”
install.packages("gap")
library(gap)
r=gcontrol2(results_short$P)
print (r$lambda)


# Using the library “GenAbel”
#install.packages("GenABEL")
library(GenABEL)
r_GenAbel <- estlambda(results_short$P, plot = FALSE, proportion = 1, method = "regression")
print(r_GenAbel$estimate)




setwd("/home/guosa/hpc/rheumatology/RA/GWAS/sample")
setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/rheumatology/RA/GWAS/sample")
RA<-read.table("RA.txt",head=T,sep="\t")
RA$FFID=unlist(lapply(RA$SID,function(x) substr(x,1,1)))
RA$LFID=unlist(lapply(RA$SID,function(x) nchar(as.character(x))))
NRA=subset(RA,FFID==3 & LFID==18)
NRA$Phen=1
NRA$BMI=NRA$Weight/(NRA$Height/100)^2
CO<-read.table("CONTROL.txt",head=T)
CO$FFID=unlist(lapply(CO$SID,function(x) substr(x,1,1)))
CO$LFID=unlist(lapply(CO$SID,function(x) nchar(as.character(x))))
NCO=subset(CO,FFID==3 & LFID==18)
NCO$Phen=0
NCO$BMI=NCO$Weight/(NCO$Height/100)^2
head(NRA)
head(NCO)
Tsex<-rbind(table(NRA$Gender),table(NCO$Gender))
Tage<-t.test(NRA$Age,NCO$Age)
Tbmi<-t.test(NRA$BMI,NCO$BMI)
par(mfrow=c(1,2))
boxplot(NCO$Age,main="Normal")
boxplot(NRA$Age,main="RA")
