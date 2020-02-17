
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

# i=match("ENSG00000213754.2",rownames(input))
# i=match("ENSG00000231246.1",rownames(input))
# i=grep("ENSG00000280054",rownames(input))
# i=grep("ENSG00000243384",rownames(input))
# i=grep("ENSG00000274370",rownames(input))
# i=grep("ENSG00000266601",rownames(input))
# i=grep("ENSG00000234883",rownames(input))

load("~/hpc/methylation/Pancancer/RNA-seq/rnaseqdata.pancancer.RData")

TCGAProjects=c("BLCA","BRCA","CESC","CHOL","COAD","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC","LUAD","LUSC","PAAD","PCPG","PRAD","READ","SARC","STAD","THCA","THYM","UCEC")
# TCGAProjects=c("BRCA")
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

idx<-which(phen$phen2==1 | phen$phen2==11)
phen<-phen[idx,]
input<-rnaseqdata[,idx]
idx<-which(phen$pid %in% paste("TCGA-",TCGAProjects,sep=""))
phen<-phen[idx,]
input<-input[,idx]
input[1:5,1:5]
input<-log(input+1,2)
input<-RawNARemove(input)
input<-RawZeroRemove(input)
Seq<-paste(phen$project_id,phen$phen2,sep="-")
rlt<-c()
coll<-c()
out1<-c()

setwd("~/hpc/meta")

preplot<-FALSE
for(i in 1:nrow(input)){
  print(i)
  mean<-tapply(as.numeric(input[i,]),Seq,function(x) mean(x,na.rm=T))
  sd<-tapply(as.numeric(input[i,]),Seq,function(x) sd(x,na.rm=T))
  num<-tapply(as.numeric(input[i,]),Seq,function(x) length(x))
  m1i=mean[seq(1,length(mean),by=2)]
  m2i=mean[seq(2,length(mean),by=2)]
  sd1i=sd[seq(1,length(mean),by=2)]
  sd2i=sd[seq(2,length(mean),by=2)]
  n1i=num[seq(1,length(mean),by=2)]
  n2i=num[seq(2,length(mean),by=2)]
  Source<-unlist(lapply(strsplit(names(m1i),"-"),function(x) x[2]))
  output<-data.frame(cbind(n1i,m1i,sd1i,n2i,m2i,sd2i))
  output$source=Source
  output<-na.omit(output)
  es<-escalc(m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i,measure="MD",data=output)
  md <- rma(es,slab=source,method = "REML", measure = "SMD",data=output)
  rlt<-rbind(rlt,c(i,md$beta,md$pval,md$ci.lb,md$ci.ub,md$I2,md$tau2))
  studlab=unlist(lapply(rownames(output),function(x) unlist(strsplit(x,"-"))[2]))
  coll<-c(coll,i)
  m<-metagen(yi,seTE=vi,data = es,comb.fixed = TRUE,comb.random = TRUE,prediction=F,sm="SMD")
  if(preplot && m$pval.random <10^-30){
  print(rownames(input)[i])
  pdf(paste(rownames(input)[i],".SMD.PANC.pdf",sep=""))
  forest(m,leftlabs =studlab,
         lab.e = "Intervention",
         pooled.totals = FALSE,
         smlab = "",studlab=studlab,
         text.random = "Overall effect",
         print.tau2 = FALSE,
         col.diamond = "blue",
         col.diamond.lines = "black",
         col.predict = "red",
         print.I2.ci = TRUE,
         digits.sd = 2,fontsize=8,xlim=c(-6,1))
  dev.off()
  }
  fixedEffect<-c(m$TE.fixed,m$lower.fixed,m$upper.fixed,m$pval.fixed)
  randomEffect<-c(m$TE.random,m$lower.random,m$upper.random,m$pval.random)
  out1<-rbind(out1,c(fixedEffect,randomEffect))
}

colnames(out1)<-c("TE.fixed","lower.fixed","upper.fixed","pval.fixed","TE.random","lower.random","upper.random","pval.random")
rownames(out1)<-row.names(input)
out1<-data.frame(out1,Symbol=ENSG2Symbol(rownames(out1)))
out1<-out1[order(out1$pval.random),]
write.table(out1,file="TCGA-SMD-DGE-Meta-Pvalue-2019.txt",sep="\t",col.names = NA,row.names = T,quote=F)
write.csv(out1,file="TCGA-SMD-DGE-Meta-Pvalue-2019.csv",quote=F)

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
out2<-c()

# i=match("ENSG00000213754.2",rownames(input))
# i=match("ENSG00000231246.1",rownames(input))
# i=grep("ENSG00000280054",rownames(input))
# i=grep("ENSG00000243384",rownames(input))
# i=grep("ENSG00000274370",rownames(input))
# i=grep("ENSG00000266601",rownames(input))
i=grep("ENSG00000234883",rownames(input))
i=grep("ENSG00000120256",rownames(input))

for(i in 1:nrow(input)){
  subpanel=FALSE
  HR<-c()
  for(TCGAProject in TCGAProjects){
    newdata<-input[,phen$project_id==paste("TCGA-",TCGAProject,sep="")]
    xphen<-phen[phen$project_id==paste("TCGA-",TCGAProject,sep=""),]
    dat<-data.frame(Rna=newdata[i,],xphen)
    thres<-mean(newdata[i,],na.rm=T)
    dat$Rna[dat$Rna<=thres]<-0
    dat$Rna[dat$Rna>thres]<-1
    hr.fit<-summary(coxph(Surv(week,censored)~Rna,dat))
    hr1=hr.fit$coefficients[1,]
    hr2=hr.fit$conf.int[1,]
    HR<-rbind(HR,c(hr1,hr2[3],hr2[4]))
    # if(subpanel==TRUE){
    # fit <- survfit(Surv(week,censored)~Rna, data = dat)
    # survp<-ggsurvplot(fit, data = dat,conf.int = F,pval = TRUE,
    #                   fun = "pct",risk.table = TRUE,size = 1,linetype = "strata",
    #                   palette = c("#E7B800","#2E9FDF"),
    #                   legend = "bottom",legend.title = rownames(input)[i],
    #                   legend.labs = c("Low-expression","High-expression"))
    # ggsave(file = paste(rownames(input)[i],"_",TCGAProject,"_KM.pdf",sep=""), survp$plot)
    # }
  }
  rownames(HR)<-TCGAProjects
  m<-metagen(HR[,1],seTE=HR[,3],comb.fixed = TRUE,comb.random = TRUE,prediction=F,sm="HR")
  # write.table(m,file=paste(rownames(input)[i],".OS.HR.PANC.txt",sep=""),sep="\t",col.names = NA,row.names = T,quote=F)
  # if(preplot && m$pval.random <10^-8){
  # pdf(paste(rownames(input)[i],".OS.HR.PANC.pdf",sep=""))
  # print(rownames(input)[i])
  # forest(m,leftlabs = rownames(HR),
  #        lab.e = "Intervention",
  #        pooled.totals = FALSE,
  #        smlab = "",studlab=rownames(HR),
  #        text.random = "Overall effect",
  #        print.tau2 = FALSE,
  #        col.diamond = "blue",
  #        col.diamond.lines = "black",
  #        col.predict = "red",
  #        print.I2.ci = TRUE,
  #        digits.sd = 2,fontsize=9)
  # dev.off()
  # }
  fixedEffect<-c(exp(m$TE.fixed),exp(m$lower.fixed),exp(m$upper.fixed),m$pval.fixed)
  randomEffect<-c(exp(m$TE.random),exp(m$lower.random),exp(m$upper.random),m$pval.random)
  out2<-rbind(out2,c(fixedEffect,randomEffect))
  # print(c(i,as.character(ENSG2Symbol(rownames(input)[i])),fixedEffect))
  print(c(i,rownames(input)[i],fixedEffect[4]))
}

colnames(out2)<-c("TE.fixed","lower.fixed","upper.fixed","pval.fixed","TE.random","lower.random","upper.random","pval.random")
rownames(out2)<-row.names(input)
out2<-data.frame(out2,ENSG2Symbol(rownames(out2)))
out3<-out2[order(out2$pval.fixed),c(1,2,3,4,9)]
out4<-na.omit(out3)
write.csv(out4,file="TCGA-HR-OS-Meta-Pvalue-2019.csv",quote=F)
write.table(out2,file="TCGA-HR-OS-Meta-Pvalue-2019.txt",sep="\t",col.names = NA,row.names = T,quote=F)
