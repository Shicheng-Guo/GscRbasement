###############################################################################################################
###   Title : Genome-wide DNA methylation analysis for atrial fibrillation patient
###   Author: Shicheng Guo, Ph.D. Email: Shihcheng.Guo@Gmail.com 
###   Section 1. Function predifinition 
###   Section 2. Data Cleaning
###   Section 3. Differential Analysis
###   Section 4. Pathway Analysis
###   Section 5. Gene Set Enrichment Analysis: enrichR in BIRC10 (https://amp.pharm.mssm.edu/Enrichr/)
###   Section 6. GEO Validation (GSE34639,GSE27895)
###   BIRC10-LC and Windows work well. CHG1 and HPC fail. 
###   11/06/2019
###############################################################################################################
BiocManager::install("ChAMP") 
BiocManager::install("doParallel") 
BiocManager::install("benchmarkme") 
BiocManager::install("enrichR") 

benchmarkme::get_ram()
detectCores()

system("cd /home/local/MFLDCLIN/guosa/hpc/methylation/Ingrid/MCaldwell-Sept27-17-HuMethEPIC/Raw_Data/idat")
setwd("/home/local/MFLDCLIN/guosa/hpc/methylation/Ingrid/MCaldwell-Sept27-17-HuMethEPIC/Raw_Data/idat")

library("ChAMP")
library("doParallel")
Dir="/home/local/MFLDCLIN/guosa/hpc/methylation/Ingrid/MCaldwell-Sept27-17-HuMethEPIC/Raw_Data/idat"
set.seed(11)
##############################################################################################################
targets <- read.metharray.sheet(Dir)
RGSet <- read.metharray.exp(targets = targets)
phenoData <- pData(RGSet)
manifest <- getManifest(RGSet)
head(getProbeInfo(manifest))
myNormalRGSet<-preprocessFunnorm(RGSet, nPCs=4, sex = NULL, bgCorr = TRUE,dyeCorr = TRUE, keepCN = TRUE, ratioConvert = TRUE,verbose = TRUE)
##############################################################################################################
##############################   Full Automatically Pipeline    ##########################################
##############################################################################################################

testDir="/home/local/MFLDCLIN/guosa/hpc/methylation/Ingrid/MCaldwell-Sept27-17-HuMethEPIC/Raw_Data/idat"
champ.process(directory = testDir)
##############################################################################################################
##############################   Separate Automatically Pipeline    ##########################################
##############################################################################################################
myLoad <- champ.load(Dir,filterBeads=TRUE,arraytype="EPIC")
# EPIC has 411 control probes
pdf("MCaldwell.AMP.EPIC.QC.pdf")
champ.QC()
dev.off()
##########################################################################
pdf("MCaldwell.AMP.EPIC.SVD.pdf")
champ.SVD(beta=myNorm,pd=myLoad$pd)
dev.off()
myCombat <- champ.runCombat(beta=myNorm,pd=myLoad$pd,batchname=c("Slide"))
##########################################################################
# don't use all the cores which will easily be killed by system
detectCores()
# seed=sample(seq(1,10000,by=1),1)
# seed=110
myNorm <- champ.norm(beta=myLoad$beta,arraytype="EPIC",cores=5,method="BMIQ")

# PureG3: 24 case vs 4 control 
myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$pureG3,compare.group=c("Case","Control"),arraytype="EPIC")
write.table(myDMP1,file=paste("AtrialFibrillation.",".24Case4Control.pureG3.myDMP.txt",sep=""),col.names = NA,row.names = T,quote=F,sep="\t")
myDMR <- champ.DMR(beta=myNorm,pheno=myLoad$pd$pureG3,method="Bumphunter",compare.group=c("Case","Control"),arraytype="EPIC",minProbes=2,cores=5,maxGap=3000)
write.table(myDMR1,file="AtrialFibrillation.24Case4Control.pureG3.myDMR.txt",col.names = NA,row.names = T,quote=F,sep="\t")
myGSEA <- champ.GSEA(beta=myNorm,DMP=myDMP[[1]], DMR=myDMR, arraytype="EPIC",adjPval=0.05, method="fisher")
myebayGSEA <- champ.ebGSEA(beta=myNorm,pheno=myLoad$pd$Sample_Group,arraytype="EPIC")
write.table(myGSEA,file="AtrialFibrillation.pureG3.myGSEA.txt",col.names = NA,row.names = T,quote=F,sep="\t")
write.table(myebayGSEA,file="AtrialFibrillation.pureG3.myebayGSEA.txt",col.names = NA,row.names = T,quote=F,sep="\t")

# Sample_Group: 24 case vs 24 control 
myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$Sample_Group,arraytype="EPIC",adjPVal = 0.01)
write.table(myDMP,file=paste("AtrialFibrillation.",".24Case24Control.myDMP.txt",sep=""),col.names = NA,row.names = T,quote=F,sep="\t")
myDMR <- champ.DMR(beta=myNorm,pheno=myLoad$pd$Sample_Group,method="Bumphunter",arraytype="EPIC",minProbes=2,cores=5,maxGap=3000)
write.table(myDMR,file="AtrialFibrillation.24Case24Control.myDMR.txt",col.names = NA,row.names = T,quote=F,sep="\t")
myGSEA <- champ.GSEA(beta=myNorm,DMP=myDMP[[1]], DMR=myDMR, arraytype="EPIC",adjPval=0.05, method="fisher")
myebayGSEA <- champ.ebGSEA(beta=myNorm,pheno=myLoad$pd$Sample_Group,arraytype="EPIC")
write.table(myGSEA$DMP,file="AtrialFibrillation.DMP.24case24control.myGSEA.txt",col.names = NA,row.names = T,quote=F,sep="\t")
write.table(myebayGSEA,file="AtrialFibrillation.24case24control.myebayGSEA.txt",col.names = NA,row.names = T,quote=F,sep="\t")
myEpiMod <- champ.EpiMod(beta=myNorm,pheno=myLoad$pd$Sample_Group)
myCNA <- champ.CNA(intensity=myLoad$intensity,pheno=myLoad$pd$Sample_Group)
myRefBase <- champ.refbase(beta=myNorm,arraytype="EPIC")
colnames(myRefBase$CorrectedBeta)=myLoad$pd$ID
rownames(myRefBase$CellFraction)=myLoad$pd$ID
write.table(myRefBase$CellFraction,file="myRefBase.CellFraction.txt",sep="\t",quote=F)
myDMPcorrect <- champ.DMP(beta = myRefBase$CorrectedBeta,pheno=myLoad$pd$Sample_Group,arraytype="EPIC",adjPVal = 1)
dim(myDMP$Case_to_Control)
P=myDMPcorrect$Case_to_Control[match(rownames(myDMP$Case_to_Control),rownames(myDMPcorrect$Case_to_Control)),]
PP=myDMP$Case_to_Control[match(rownames(P[which(p.adjust(P$P.Value,"bonferroni")<0.05),]),rownames(myDMP$Case_to_Control)),]
write.table(PP,file="Table_1.txt",col.names = NA,row.names = T,quote=F,sep="\t")
myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$Sample_Group,arraytype="EPIC",adjPVal = 1)
ManhattanmyDMP(myDMP$Case_to_Control)

# Young_Old: 24 case vs 24 control 
myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$Young_Old,arraytype="EPIC",adjPVal = 0.05)
write.table(myDMP,file="AtrialFibrillation.YoungOld.myDMP.txt",col.names = NA,row.names = T,quote=F,sep="\t")
myDMR <- champ.DMR(beta=myNorm,pheno=myLoad$pd$Young_Old,method="Bumphunter",compare.group=c("Old","Young"),arraytype="EPIC",minProbes=2,cores=6,maxGap=300)
write.table(myDMR,file="AtrialFibrillation.YoungOld.myDMR.txt",col.names = NA,row.names = T,quote=F,sep="\t")

myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$Young_Old,arraytype="EPIC",adjPVal = 1)
ManhattanmyDMP(myDMP$Old_to_Young)
myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$cmG1,arraytype="EPIC",adjPVal =1)
ManhattanmyDMP(myDMP$Old_Case_to_Old_Control)
ManhattanmyDMP(myDMP$Young_Case_to_Young_Control)

myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$cmG1,arraytype="EPIC",adjPVal =0.05)
write.table(myDMP$Old_Case_to_Old_Control,file="AtrialFibrillation.Old_Case_to_Old_Control.myDMR.txt",col.names = NA,row.names = T,quote=F,sep="\t")
write.table(myDMP$Young_Case_to_Young_Control,file="AtrialFibrillation.Young_Case_to_Young_Control.myDMR.txt",col.names = NA,row.names = T,quote=F,sep="\t")
myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$cmG2,arraytype="EPIC",adjPVal =0.05)
write.table(myDMP$Case_Before_to_Control_Before,file="AtrialFibrillation.Case_Before_to_Control_Before.myDMR.txt",col.names = NA,row.names = T,quote=F,sep="\t")
myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$cmG2,arraytype="EPIC",adjPVal =1)
ManhattanmyDMP(myDMP$Case_Before_to_Control_Before)


x<-head(myLoad$beta)
library('caret')
corr<-cor(t(myNorm))
setwd("~/hpc/methylation/Ingrid/MCaldwell-Sept27-17-HuMethEPIC/Raw_Data/idat")
save(myNorm,file="myNorm.RData")

corr<-c()
for(i in seq(1000,200000,5000)){
x<-sample(1:nrow(myNorm),i)
corr<-c(corr,median(cor(t(myNorm[x,]))))
print(i)
}

map<-read.csv("~/hpc/db/hg19/GPL21145_MethylationEPIC_15073387_v-1-0.csv",skip=7,head=T)
newmap<-map[,c(12,13,1,7,15,18,19)]
newmap$CHR=paste("chr",newmap$CHR,sep="")
newmap$START=newmap$MAPINFO-1
newmap$END=newmap$MAPINFO
bed<-newmap[,c(1,8,9,3)]
write.table(bed,file="MethylationEPIC.hg19.bed",sep="\t",quote=F,col.names = F,row.names = F)

map<-read.table("~/hpc/db/hg19/MethylationEPIC.hg19.sort.bed",head=F)

input<-t(myNorm[na.omit(match(subset(map,map[1]=="chr22")[,4],rownames(myNorm))),])
corr<-cor(input)
pdf("hist.pdf")
hist(corr,main="Histogram of correlation between variables",col="blue")
dev.off()

hc<-findCorrelation(corr,cutoff=0.4)
length(hc)

# estimate the best cor cutoff threshold for findCorrelation() for EPIC dataset (Atrial Fibrillation)
input<-t(myNorm[na.omit(match(subset(map,map[1]=="chr22")[,4],rownames(myNorm))),])
R<-c()
P<-c()
for(i in 1:ncol(input)){
  print(i)
  for(j in (i+1):ncol(input)){
    x<-input[,i]
    y<-input[,j]
    fit<-cor.test(x,y)
    r=abs(fit$estimate)
    p=fit$p.value
    R<-c(R,r)
    P<-c(P,p)
  }
}
rlt<-list()
rlt$R<-R
rlt$P<-P
save(rlt,file="P-R-AF-2019-Shicheng-Guo.RData")
pdf("P-R-AF-2019-Shicheng-Guo.pdf")
plot(x=R,y=-log(P,10),pch=16)
dev.off()


## DMR analysis
myDMR <- champ.DMR(beta=myNorm,pheno=myLoad$pd$cmG1,method="Bumphunter",compare.group=c("Old_Case","Old_Control"),arraytype="EPIC",minProbes=2,cores=5,maxGap=300)
write.table(myDMR,file="AtrialFibrillation.Old_Case_to_Old_Control.myDMR.txt",col.names = NA,row.names = T,quote=F,sep="\t")
myDMR <- champ.DMR(beta=myNorm,pheno=myLoad$pd$cmG1,method="Bumphunter",compare.group=c("Old_Case","Young_Control"),arraytype="EPIC",minProbes=2,cores=5,maxGap=300)
write.table(myDMR,file="AtrialFibrillation.Old_Case_to_Young_Control.myDMR.txt",col.names = NA,row.names = T,quote=F,sep="\t")
myDMR <- champ.DMR(beta=myNorm,pheno=myLoad$pd$cmG1,method="Bumphunter",compare.group=c("Young_Case","Young_Control"),arraytype="EPIC",minProbes=2,cores=5,maxGap=300)
write.table(myDMR,file="AtrialFibrillation.Young_Case_to_Young_Control.myDMR.txt",col.names = NA,row.names = T,quote=F,sep="\t")


myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$Pmr_Enrollment,arraytype="EPIC",adjPVal =0.05)
ManhattanmyDMP(myDMP$After_to_Before)
write.table(myDMP,file="AtrialFibrillation.Enrollment.myDMP.txt",col.names = NA,row.names = T,quote=F,sep="\t")
myDMR <- champ.DMR(beta=myNorm,pheno=myLoad$pd$Pmr_Enrollment,method="Bumphunter",compare.group=c("Case","Control"),arraytype="EPIC",minProbes=2,cores=6,maxGap=300)
write.table(myDMR,file="AtrialFibrillation.Enrollment.myDMR.txt",col.names = NA,row.names = T,quote=F,sep="\t")

myGSEA <- champ.GSEA(beta=myNorm,DMP=myDMP[[1]], DMR=myDMR, arraytype="EPIC",adjPval=0.05, method="fisher")
myebayGSEA <- champ.ebGSEA(beta=myNorm,pheno=myLoad$pd$Sample_Group,arraytype="EPIC")
write.table(myebayGSEA,file="AtrialFibrillation.myebayGSEA.txt",col.names = NA,row.names = T,quote=F,sep="\t")

myEpiMod <- champ.EpiMod(beta=myNorm,pheno=myLoad$pd$Sample_Group)
myCNA <- champ.CNA(intensity=myLoad$intensity,pheno=myLoad$pd$Sample_Group)
myRefBase1 <- champ.refbase(beta=myNorm,arraytype="EPIC")

myEpiMod <- champ.EpiMod(beta=myNorm,pheno=myLoad$pd$Sample_Group)
myCNA <- champ.CNA(intensity=myLoad$intensity,pheno=myLoad$pd$Sample_Group)
myRefBase1 <- champ.refbase(beta=myNorm,arraytype="EPIC")

myBlock <- champ.Block(beta=myNorm,pheno=myLoad$pd$Sample_Group,arraytype="EPIC",cores=6)
write.table(myBlock$Block,file="AtrialFibrillation.CaseControl.myblock.txt",col.names = NA,row.names = T,quote=F,sep="\t")
myBlock <- champ.Block(beta=myNorm,pheno=myLoad$pd$Sample_Group,arraytype="EPIC",cores=6)
write.table(myBlock$Block,file="AtrialFibrillation.YoungOld.myblock.txt",col.names = NA,row.names = T,quote=F,sep="\t")
myBlock <- champ.Block(beta=myNorm,pheno=myLoad$pd$Sample_Group,arraytype="EPIC",cores=6)
write.table(myBlock$Block,file="AtrialFibrillation.Enrollment.myblock.txt",col.names = NA,row.names = T,quote=F,sep="\t")

ManhattanmyDMP<-function(myDMP){
  library("qqman")
  library("Haplin")
  SNP=rownames(myDMP)
  CHR=myDMP$CHR
  if(length(grep("X",CHR))>0){
    CHR<-sapply(CHR,function(x) gsub(pattern = "X",replacement = "23",x))
    CHR<-sapply(CHR,function(x) gsub(pattern = "Y",replacement = "24",x))
  }
  CHR<-as.numeric(CHR)
  BP=myDMP$MAPINFO
  P=myDMP$P.Value
  manhattaninput=na.omit(data.frame(SNP,CHR,BP,P))
  max<-max(2-log(manhattaninput$P,10))
  genomewideline=0.05/nrow(manhattaninput)
  seed=sample(1:10000,1)
  pdf(paste("manhattan.",seed,".pdf",sep=""))
  manhattan(manhattaninput,col = c("blue4", "orange3"),ylim = c(0,10),lwd=2, suggestiveline=F,genomewideline=FALSE)
  dev.off()
  pdf(paste("qqplot",seed,".pdf",sep=""))
  pQQ(P, nlabs =length(P), conf = 0.95)
  dev.off()
}

 
