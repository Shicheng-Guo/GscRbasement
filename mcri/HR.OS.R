source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/GscTools.R")
library("meta")
library("metafor")
library("survival")
library("survminer")

Symbol2ENSG<-function(Symbol){
  db<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/ENSG.ENST.ENSP.Symbol.hg19.bed",sep="\t")
  ENSG<-as.character(db[match(Symbol,db$V4),8])
  ENSG<-na.omit(data.frame(Symbol,ENSG))
  return(ENSG)
}
ENSG2Symbol<-function(ENSG){
  db<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/ENSG.ENST.ENSP.Symbol.hg19.bed",sep="\t")
  ENSG<-unlist(lapply(strsplit(ENSG,split="[.]"),function(x) x[1]))
  Symbol<-db[match(as.character(ENSG),db$V8),4]
  return(Symbol)
}

ensg2bed<-function(ENSG){
  db<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/hg19/ENSG.ENST.hg19.txt",as.is=T,head=F)
  ENSG<-unlist(lapply(strsplit(ENSG,split="[.]"),function(x) x[1]))
  bed<-unique(db[db$V5 %in% as.character(ENSG),c(1,2,3,5)])
  return(bed)
}

load("~/hpc/methylation/Pancancer/RNA-seq/rnaseqdata.pancancer.RData")

TCGAProjects=c("BLCA","BRCA","CESC","CHOL","COAD","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC","LUAD","LUSC","PAAD","PCPG","PRAD","READ","SARC","STAD","THCA","THYM","UCEC")
panc<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/PANC/master/extdata/panc.txt",head=T)
phen1=read.table("https://raw.githubusercontent.com/Shicheng-Guo/PANC/master/extdata/TCGA-clinical-11093.tsv",header = T,sep="\t")
phen2=read.table("https://raw.githubusercontent.com/Shicheng-Guo/PANC/master/extdata/File_metadata2.txt",header = T,sep="\t")
head(phen1)
head(phen2)
colnames(rnaseqdata)<-unlist(lapply(strsplit(colnames(rnaseqdata),"/"),function(x) x[2]))
phen<-data.frame(phen2,phen1[match(phen2$cases.0.case_id,phen1$case_id),])
phen$file_name=gsub(".gz","",phen$file_name)
# prepare phenotype information
phen<-phen[match(colnames(rnaseqdata),phen$file_name),]
phen$phen4<-id2phen4(phen$cases.0.samples.0.submitter_id)
phen$phen3<-id2phen3(phen$cases.0.samples.0.submitter_id)
phen$phen2<-id2bin(phen$cases.0.samples.0.submitter_id)
phen$pid<-phen$project_id
head(phen)

OS<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/OverallSurvivalTime.txt",head=T,sep="\t")
# match survival information
idx<-which(c(phen$phen2==1))
phen<-phen[idx,]
input<-rnaseqdata[,idx]
input[1:5,1:5]
idx<-na.omit(match(OS$submitter_id,phen$phen3))
input<-input[,idx]
phen<-phen[idx,]
phen<-data.frame(phen,OS[match(phen$phen3,OS$submitter_id),])
phen$censored<-as.numeric(!phen$censored)
phen$week=phen$time/7

i=grep("ENSG00000231246",rownames(input))          # PANC246
i=grep("ENSG00000213754",rownames(input))          # PANC754
i=grep("ENSG00000181896",rownames(input))          # ZNF101
i=grep("ENSG00000131849",rownames(input))          # ZNF132
i=grep(Symbol2ENSG("NCOA4")[1,2],rownames(input))  # NCOA4
i=grep(Symbol2ENSG("SLC7A11")[1,2],rownames(input))# SLC7A11
i=grep(Symbol2ENSG("SAT2")[1,2],rownames(input))   # SAT2
i=grep(Symbol2ENSG("GSS")[1,2],rownames(input))    # GSS
i=grep("ENSG00000120256",rownames(input))          # LRP11


xii<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/ferroptosis/master/codependency.txt",as.is=T)[,1]
xii<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/ferroptosis/master/ferroptosis.genelist.txt",sep="\t",as.is=T)[,2]
xii<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/ferroptosis/master/SLC7A11.codependency.txt",as.is=T)[,1]
xii<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/breast/master/extdata/brca.tcga.target.hg19.bed",sep="\t",as.is=T)[,4]
xii<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/cholangiocarcinoma/master/cholangiocarcinoma.hg19.bed",sep="\t",as.is=T)[,4]
xii<-rownames(input)

setwd("/home/guosa/hpc/meta/OS")

out2<-c()
for(i in 1:nrow(input)){
HR<-c()
Theshold<-c()
for(TCGAProject in TCGAProjects){
  newdata<-input[,phen$project_id==paste("TCGA-",TCGAProject,sep="")]
  xphen<-phen[phen$project_id==paste("TCGA-",TCGAProject,sep=""),]
  dat<-data.frame(Rna=log(newdata[i,]+1,2),xphen)
  thres<-mean(dat[,1],na.rm=T)
  Theshold<-c(Theshold,thres)
  #dat$Rna[dat$Rna<=thres]<-0
  #dat$Rna[dat$Rna>thres]<-1
  hr.fit<-summary(coxph(Surv(week,censored)~Rna,dat))
  hr1=hr.fit$coefficients[1,]
  hr2=hr.fit$conf.int[1,]
  HR<-rbind(HR,c(hr1,hr2[3],hr2[4]))
  
  fit <- survfit(Surv(week,censored)~Rna, data = dat)
  # survp<-ggsurvplot(fit, data = dat,conf.int = F,pval = TRUE,
  #                 fun = "pct",risk.table = TRUE,size = 1,linetype = "strata",
  #                 palette = c("#2E9FDF","#E7B800"),
  #                 legend = "bottom",legend.title = rownames(input)[i],
  #                 legend.labs = c("Low-expression","High-expression"))
  # ggsave(file = paste(ENSG2Symbol(rownames(input)[i]),"-",rownames(input)[i],"_",TCGAProject,"_KM.pdf",sep=""), survp$plot)
}

print(i)
rownames(HR)<-TCGAProjects
m<-metagen(HR[,1],seTE=HR[,3],comb.fixed = TRUE,comb.random = TRUE,prediction=F,sm="HR")
fixedEffect<-c(exp(m$TE.fixed),exp(m$lower.fixed),exp(m$upper.fixed),m$pval.fixed)
randomEffect<-c(exp(m$TE.random),exp(m$lower.random),exp(m$upper.random),m$pval.random)
out2<-rbind(out2,c(fixedEffect,randomEffect))

pdf(paste(ENSG2Symbol(rownames(input)[i]),"-",rownames(input)[i],".OS.HR.PANC.pdf",sep=""))
print(rownames(input)[i])
forest(m,leftlabs = rownames(HR),
       lab.e = "Intervention",
       pooled.totals = FALSE,
       smlab = "",studlab=rownames(HR),
       text.random = "Overall effect",
       print.tau2 = FALSE,
       col.diamond = "blue",
       col.diamond.lines = "black",
       col.predict = "red",
       print.I2.ci = TRUE,
       digits.sd = 2,fontsize=9,xlim=c(0.2,3))
dev.off()
write.table(HR,file=paste(ENSG2Symbol(rownames(input)[i]),"-",rownames(input)[i],".OS.HR.EACH.txt",sep=""),sep="\t",quote=F,col.names=NA,row.names=T)
}
colnames(out2)<-c("TE.fixed","lower.fixed","upper.fixed","pval.fixed","TE.random","lower.random","upper.random","pval.random")
# rownames(out2)<-as.character(ENSG2Symbol(as.character(Symbol2ENSG(xii)[,2])))
rownames(out2)<-rownames(input)
write.csv(out2,file=paste("Cancer.methylation.targets","OS.HR.csv",sep=""),quote=F)
