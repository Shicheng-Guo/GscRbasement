data<-read.table("ROI.RA3000.model",head=T,as.is=T)

ctabe<-function(input){
input<-as.numeric(unlist(strsplit(as.character(input),"/")))
geno.n<-input[c(3,2,1)]
geno.p<-percent(input[c(3,2,1)]/sum(input[1:3]))
dom.n<-input[1]+input[2]
dom.p<-percent((input[1]+input[2])/sum(input[1:3]))
rec.n<-input[1]
rec.p<-percent(input[1]/sum(input[1:3]))
risk.n<-2*input[1]+input[2]
risk.p<-percent((2*input[1]+input[2])/(2*sum(input[1:3])))
n<-c(geno.n,rec.n,dom.n,risk.n)
p<-c(geno.p,rec.p,dom.p,risk.p)
return(c(paste(n,"(",p,")",sep="")))
}

library("epitools")
outable<-c()
for(i in unique(data$SNP)){
temp<-subset(data,SNP==i)
mia<-temp$A1[1]
maa<-temp$A2[1]
AA<-paste(mia,mia,sep="")
AC<-paste(mia,maa,sep="")
CC<-paste(maa,maa,sep="")
rec<-paste(AA," vs ",AC,"+",CC,sep="")
dom<-paste(AA,"+",AC," vs ",CC,sep="")
A<-paste(mia,"allele",sep=" ")
TYPE<-c(CC,AC,AA,dom,rec,A)

AFF=ctabe(temp$AFF)
UNAFF=ctabe(temp$UNAFF)
data.frame(AFF,UNAFF)

REF1<-as.numeric(unlist(lapply(strsplit(as.character(AFF),split="[(]"),function(x) x[1]))[1:6])
REF2<-as.numeric(unlist(lapply(strsplit(as.character(UNAFF),split="[(]"),function(x) x[1]))[1:6])

REF<-data.frame(REF1,REF2)
OR<-"1 (reference)"
P<-" "
for(j in 2:3){
xx<-as.matrix(REF[c(1,j),])
or<-round((xx[2,1]/xx[1,1])/(xx[2,2]/xx[1,2]),2)
ci95<-round(oddsratio(xx, y = NULL,method = c("fisher"), rev="col",conf.level = 0.95,correction = T,verbose = FALSE)$measure[2,2:3],2)
or<-paste(or,"(",ci95[1],"-",ci95[2],")",sep="")
p<-fisher.test(xx)$p.value
OR<-c(OR,or)
P<-c(P,p)
}
for(j in 4:5){
xx<-matrix(c(sum(REF[1:3,1])-REF[j,1],REF[j,1],sum(REF[1:3,2])-REF[j,2],REF[j,2]),2,2,byrow=F)
or<-round((xx[2,1]/xx[1,1])/(xx[2,2]/xx[1,2]),2)
ci95<-round(oddsratio(xx, y = NULL,method = c("fisher"), rev="col",conf.level = 0.95,correction = T,verbose = FALSE)$measure[2,2:3],2)
or<-paste(or,"(",ci95[1],"-",ci95[2],")",sep="")
p<-fisher.test(xx)$p.value
OR<-c(OR,or)
P<-c(P,p)
}
