# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Fri Feb 1 15:28:27 EST 2019

################################################################
#   Differential expression analysis with limma
install.packages("purrr")
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE76809", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL6480", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXX1111111111111111111100000000000000000000000000",
               "00000000000000000000000000000000000001100000000000",
               "00000000000000000000000000000000000000000000000000",
               "0000000000000000000XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXX000000000011111")
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)

####################################
library(stats)
library(edgeR)
library("RCurl")
library("qualityTools")
library(devtools)
library(Biobase)
library(preprocessCore)

input<-exprs(gset)
input[1:5,1:5]
newinput = normalize.quantiles(as.matrix(input))
colnames(newinput)<-colnames(input)
rownames(newinput)<-rownames(input)
hkurl <- getURL("https://raw.githubusercontent.com/CNAID/mRNA/master/Human_Housekeeping_Gene_List.txt")
hk <- read.table(text = hkurl)
GPL6480<-fData(gset)
PRID<-GPL6480[na.omit(match(hk[,2],GPL6480$GenBank.Accession)),]$ID
input.control<-newinput[match(PRID,rownames(newinput)),]
pca<-prcomp(x=na.omit(input.control))
confounders<-pca$rotation[,1:1]
levels(fl)<-c(1,0)
table(fl)
rlt<-matrix(rep(NA,3*nrow(input)),ncol=3)
for(i in 1:nrow(input)){
print(i)
NewInput<-na.omit(data.frame(phen=as.numeric(as.character(fl)),x=newinput[i,],confounders))
head(NewInput)
if(nrow(NewInput)>2 && var(NewInput$x)>0){
try(glmrlt<-summary(glm(phen~.,NewInput,family =gaussian)))
rlt[i,1:3]<-c(i,glmrlt$coefficients[2,3:4])
}
}
colnames(rlt)<-c("iid","t-value","pvals")
output<-cbind(rlt,GPL6480)

head(output)
head(hk)
hk.ssc<-output[na.omit(match(hk[,2],output$GenBank.Accession)),]
par(mfrow=c(2,2))
qqplot(na.omit(output$pvals))
qqplot(na.omit(hk.ssc$pvals))

ironurl <- getURL("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/FGF6/Iron.Related.Gene.txt")
iron <- read.table(text = ironurl)
IronID<-na.omit(GPL6480[match(iron[,1],GPL6480$Gene.symbol),]$ID)
ironpvals<-output[match(IronID,output$ID),]
qqplot(na.omit(ironpvals$pvals))

write.table(output,file="GEO-SSc-Vs-Normal-GPL6480.Pvalue.txt",sep="\t",quote=F,row.names = T,col.names = NA)

####################################
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

library(QQperm)
lambda <-estlambda2(example.Ps$observed,example.Ps$perm, plot = TRUE, adjust.xy = TRUE)

head(example.Ps)

qqPlot(x, "logistic")   

################################################################
#   Boxplot for selected GEO samples
library(Biobase)
library(GEOquery)

# load series and platform data from GEO

gset <- getGEO("GSE76809", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL6480", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# group names for all samples in a series
gsms <- paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXX1111111111111111111100000000000000000000000000",
               "00000000000000000000000000000000000001100000000000",
               "00000000000000000000000000000000000000000000000000",
               "0000000000000000000XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXX000000000011111")
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
sml <- paste("G", sml, sep="")  

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# order samples by group
ex <- exprs(gset)[ , order(sml)]
sml <- sml[order(sml)]
fl <- as.factor(sml)
labels <- c("SSc_Skin","Normal_Skin")

# set parameters and draw the plot
palette(c("#dfeaf4","#f4dfdf", "#AABBCC"))
dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE76809", '/', annotation(gset), " selected samples", sep ='')
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, fill=palette(), bty="n")
