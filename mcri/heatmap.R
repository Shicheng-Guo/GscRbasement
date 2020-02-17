df <- scale(mtcars)
heatmap(x, scale = "row")
heatmap(df, scale = "none")
col<- colorRampPalette(c("red", "white", "blue"))(256)
library("RColorBrewer")
col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
library("RColorBrewer")
col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
heatmap(df, scale = "none", col =  col, 
        RowSideColors = rep(c("blue", "pink"), each = 16),
        ColSideColors = c(rep("purple", 5), rep("orange", 6)))
# install.packages("gplots")
library("gplots")
heatmap.2(df, scale = "none", col = bluered(100), 
          trace = "none", density.info = "none")

# Pretty heat maps: pheatmap()
install.packages("pheatmap")
library("pheatmap")
pheatmap(df, cutree_rows = 4)

# Interactive heat maps: d3heatmap()
library("d3heatmap")
d3heatmap(scale(mtcars), colors = "RdYlBu",k_row = 4, k_col = 2)


library(dendextend)
Rowv  <- mtcars %>% scale %>% dist %>% hclust %>% as.dendrogram %>% set("branches_k_color", k = 3) %>% set("branches_lwd", 1.2) %>%ladderize
Colv  <- mtcars %>% scale %>% t %>% dist %>% hclust %>% as.dendrogram %>%set("branches_k_color", k = 2, value = c("orange", "blue")) %>%set("branches_lwd", 1.2) %>%ladderize

heatmap(scale(mtcars), Rowv = Rowv, Colv = Colv,scale = "none")

library(gplots)
heatmap.2(scale(mtcars), scale = "none", col = bluered(100), Rowv = Rowv, Colv = Colv,trace = "none", density.info = "none")

library("d3heatmap")
d3heatmap(scale(mtcars), colors = "RdBu",Rowv = Rowv, Colv = Colv)



BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
df <- scale(mtcars)
Heatmap(df, 
        name = "mtcars", #title of legend
        column_title = "Variables", row_title = "Samples",
        row_names_gp = gpar(fontsize = 7) # Text size for row names
)

library(circlize)
mycols <- colorRamp2(breaks = c(-2, 0, 2),colors = c("green", "white", "red"))
Heatmap(df, name = "mtcars", col = mycols)

library("circlize")
library("RColorBrewer")
Heatmap(df, name = "mtcars",col = colorRamp2(c(-2, 0, 2), brewer.pal(n=3, name="RdBu")))


library("dendextend")
row_dend = hclust(dist(df)) # row clustering
col_dend = hclust(dist(t(df))) # column clustering
Heatmap(df, name = "mtcars", 
        row_names_gp = gpar(fontsize = 6.5),
        cluster_rows = color_branches(row_dend, k = 4),
        cluster_columns = color_branches(col_dend, k = 2))

set.seed(2)
Heatmap(df, name = "mtcars", k = 2)

Heatmap(df, name = "mtcars", split = mtcars$cyl,row_names_gp = gpar(fontsize = 7))

# Split by combining multiple variables
Heatmap(df, name ="mtcars", split = data.frame(cyl = mtcars$cyl, am = mtcars$am),row_names_gp = gpar(fontsize = 7))

Heatmap(df, name ="mtcars", col = mycol,km = 2, split =  mtcars$cyl)

# install.packages("cluster")
library("cluster")
set.seed(2)
df <- scale(mtcars)
pa = pam(df, k = 3)
mycol <- colorRamp2(breaks = c(-2, 0, 2),colors = c("green", "white", "red"))

Heatmap(df, name = "mtcars", col = mycol,split = paste0("pam", pa$clustering))
HeatmapAnnotation(df, name, col, show_legend)
df <- t(df)
? HeatmapAnnotation
# Annotation data frame
df <- data.frame(cyl = mtcars$cyl, am = mtcars$am,mpg = mtcars$mpg)
# Define colors for each levels of qualitative variables
# Define gradient color for continuous variable (mpg)
col = list(cyl = c("4" = "green", "6" = "gray", "8" = "darkred"),
           am = c("0" = "yellow", "1" = "orange"),
           mpg = circlize::colorRamp2(c(17, 25), c("lightblue", "purple"))
           )
# Create the heatmap annotation
ha <- HeatmapAnnotation(df,col)
# Combine the heatmap and the annotation
Heatmap(df, name = "mtcars",top_annotation = ha)



# Define some graphics to display the distribution of columns
.hist = anno_histogram(df, gp = gpar(fill = "lightblue"))
.density = anno_density(df, type = "line", gp = gpar(col = "blue"))
ha_mix_top = HeatmapAnnotation(hist = .hist, density = .density)
# Define some graphics to display the distribution of rows
.violin = anno_density(df, type = "violin", 
                       gp = gpar(fill = "lightblue"), which = "row")
.boxplot = anno_boxplot(df, which = "row")
ha_mix_right = HeatmapAnnotation(violin = .violin, bxplt = .boxplot,
                                 which = "row", width = unit(4, "cm"))
# Combine annotation with heatmap
Heatmap(df, name = "mtcars", 
        column_names_gp = gpar(fontsize = 8),
        top_annotation = ha_mix_top, 
        top_annotation_height = unit(3.8, "cm")) + ha_mix_right
# Heatmap 1
ht1 = Heatmap(df, name = "ht1", km = 2,
              column_names_gp = gpar(fontsize = 9))
# Heatmap 2
ht2 = Heatmap(df, name = "ht2", 
              col = circlize::colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
              column_names_gp = gpar(fontsize = 9))
# Combine the two heatmaps
ht1 + ht2

install.packages("goeveg")
install.packages("yaml")
library("RColorBrewer")
library("gplots")
library("pheatmap")
library("d3heatmap")
library("dendextend")
library("ComplexHeatmap")
library("circlize")
library("goeveg")

percent<-function(x, digits = 2, format = "f", ...){
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}

set.seed(123)
mat = matrix(rnorm(2000), nrow=200,ncol=10)
head(mat)
rownames(mat) = paste("R", 1:200,sep="")
colnames(mat) = paste("C", 1:10,sep="")
column_ha = HeatmapAnnotation(foo1 = runif(10), bar1 = anno_barplot(runif(10)))
row_ha = rowAnnotation(foo2 = runif(200), bar2 = anno_barplot(runif(200)))
Heatmap(mat, name = "mat", top_annotation = column_ha, right_annotation = row_ha)


setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/RnaseqBacterial/extdata/rnaseq/Rockhopper_Results")

perl -p -i -e "s/\'/-/g" gene.tab.matrix.v8.txt
perl -p -i -e "s/\(/-/g" gene.tab.matrix.v8.txt
perl -p -i -e "s/\)/-/g" gene.tab.matrix.v8.txt
perl -p -i -e "s/\,/-/g" gene.tab.matrix.v8.txt
perl -p -i -e "s/\:/-/g" gene.tab.matrix.v8.txt
perl -p -i -e "s/\;/-/g" gene.tab.matrix.v8.txt

data<-read.table("gene.tab.matrix.v8.txt",head=T,row.names=1,sep="\t",as.is=T,check.names = F)
map<-read.table("gene.tab.matrix.v8.txt",head=T,row.names=1,sep="\t",as.is=T,check.names = F)
head(map)

data<-data[order(apply(data[,1:24],1,function(x) sum(x>0)),decreasing = T),]
data<-data[,c(order(as.numeric(colnames(data[1:24]))),25,26)]

data$pc<-percent(apply(data[,1:24],1,function(x) sum(x>0))/24)
data$mean<-round(apply(data[,1:24],1,function(x) mean(x)),0)
data$sd<-round(apply(data[,1:24],1,function(x) sd(x)),0)
data$cv<-round(apply(data[,1:24],1,function(x) cv(x)),2)

head(data)
write.csv(data,file="gene.tab.matrix.v9.csv",quote=F,row.names = T)

GOI1<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/RnaseqBacterial/master/extdata/interestlist1.txt")
GOI2<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/RnaseqBacterial/master/extdata/interestlist2.txt")

x1<-tolower(GOI1[,1])
x2<-tolower(GOI1[,1])
x3<-tolower(data$N1)

x1[which( ! x1 %in% x3)]
x2[which( ! x2 %in% x3)]

out1<-na.omit(data[match(tolower(GOI1[,1]),tolower(data$N1)),])
out2<-na.omit(data[match(tolower(GOI2[,1]),tolower(data$N1)),])
out1
out2
write.csv(out1,file="gene.tab.matrix.table3.csv",quote=F)
write.csv(out2,file="gene.tab.matrix.table4.csv",quote=F)


sam<-read.csv("https://raw.githubusercontent.com/Shicheng-Guo/RnaseqBacterial/master/extdata/saminfo.csv")
sam<-sam[na.omit(match(colnames(data),sam$SID)),]
write.csv(sam,file="duRNAseq_sample_information.csv",quote=F,row.names=F)

data<-data[match(unique(rbind(out1,out2))$N1,data$N1),]
rownames(data)<-data$N1
head(data)

data<-data[,-match(c(8,45),colnames(data))]
data<-data[,1:22]
anno<-sam[,c(4,5,6,7,8,12)]
dim(anno)
table(anno$ST)

#######################
## heatmap analysis
######################
library("gplots")
library("circlize")
library("ComplexHeatmap")
set.seed(123)
data<-read.csv(file="https://raw.githubusercontent.com/Shicheng-Guo/RnaseqBacterial/master/extdata/RPKM.dual-RNA-seq.csv",row.names = 1,as.is = T)
anno<-read.csv(file="https://raw.githubusercontent.com/Shicheng-Guo/RnaseqBacterial/master/extdata/Saminfo.RPKM.dual-RNA-seq.csv",row.names = 1,as.is = T)
anno<-read.csv(file="https://raw.githubusercontent.com/Shicheng-Guo/RnaseqBacterial/master/extdata/saminfo.csv",row.names = 1,as.is = T)
anno
annodf<-data.frame(cid=anno$CID,gender=anno$Gender,age=anno$Age,severity=anno$Severity)
column_ha = HeatmapAnnotation(CID=anno$CID,Gender=anno$Gender,Age=anno$Age,Severity=anno$Severity,ST=anno$ST)
Heatmap(scale(data.matrix(data),center = T,scale=F),name="RPKM",top_annotation = column_ha)
column_ha = HeatmapAnnotation(CID=anno$CID,Gender=anno$Gender,Age=anno$Age,Severity=anno$Severity,ST=anno$ST)
Heatmap(data.matrix(data),name="RPKM",top_annotation = column_ha)
#######################
######################

library("arm")
RowNARemove<-function(data,missratio=0.3){
  threshold<-(missratio)*dim(data)[2]
  NaRaw<-which(apply(data,1,function(x) sum(is.na(x))>threshold))
  zero<-which(apply(data,1,function(x) all(x==0))==T)
  NaRAW<-c(NaRaw,zero)
  if(length(NaRAW)>0){
    dat<-data[-NaRAW,]
  }else{
    dat<-data;
  }
  dat
} 

data<-RowNARemove(data)
dim(data)
OR<-exp(apply(data,1,function(x) summary(bayesglm(as.numeric(as.factor(anno$Severity))~as.numeric(x)))$coefficients[2,1]))
P<-apply(data,1,function(x) summary(bayesglm(as.numeric(as.factor(anno$Severity))~as.numeric(x)))$coefficients[2,4])
Risk<-data.frame(OR,P)
head(Risk)
input<-Risk
dim(input)
library("CMplot")
Position=map[match(rownames(input),rownames(map)),]$POS
cminput<-na.omit(data.frame(SNP=rownames(input),Chromosome=rep(1,nrow(input)),Position=Position,trait1=input$P))
RawPos=cminput$Position
bin=round(max(cminput$Position)/22)
bin
cminput$Chromosome=round(cminput$Position/bin)
cminput$Position=cminput$Position-(cminput$Chromosome-1)*bin
CMplot(cminput,plot.type="b",LOG10=TRUE,file="jpg",memo="Severity",threshold.lty=c(1,2),signal.col=c("red","green"),dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6)
cminput$N1=map[match(cminput$SNP,rownames(map)),]$N1
cminput$N2=map[match(cminput$SNP,rownames(map)),]$N2
cminput$Rawposition=RawPos
cminput<-cminput[order(cminput$trait1,decreasing = F),]
write.table(cminput,file="Severity.manhattan.input.txt",sep="\t",quote=F,col.names =T,row.names = F)
write.csv(Risk,file="Severity.Risk.OR.Pvalue.csv",quote=F)


OR<-exp(apply(data,1,function(x) summary(bayesglm(as.numeric(as.factor(anno$ST))~as.numeric(x)))$coefficients[2,1]))
P<-apply(data,1,function(x) summary(bayesglm(as.numeric(as.factor(anno$ST))~as.numeric(x)))$coefficients[2,4])
Risk<-data.frame(OR,P)
head(Risk)
input<-Risk
library("CMplot")
Position=map[match(rownames(input),rownames(map)),]$POS
cminput<-na.omit(data.frame(SNP=rownames(input),Chromosome=rep(1,nrow(input)),Position=Position,trait1=input$P))
RawPos=cminput$Position
bin=round(max(cminput$Position)/22)
bin
cminput$Chromosome=round(cminput$Position/bin)
cminput$Position=cminput$Position-(cminput$Chromosome-1)*bin
CMplot(cminput,plot.type="b",LOG10=TRUE,file="jpg",memo="PubSwab",threshold.lty=c(1,2),signal.col=c("red","green"),dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6)
cminput$N1=map[match(cminput$SNP,rownames(map)),]$N1
cminput$N2=map[match(cminput$SNP,rownames(map)),]$N2
cminput$Rawposition=RawPos
cminput<-cminput[order(cminput$trait1,decreasing = F),]
write.table(cminput,file="PubSwab.manhattan.input.txt",sep="\t",quote=F,col.names =T,row.names = F)
write.csv(Risk,file="Pubs.Swab.Risk.OR.Pvalue.csv",quote=F)

getwd()
OR<-exp(apply(data,1,function(x) summary(bayesglm(as.numeric(as.factor(anno$Gender))~as.numeric(x)))$coefficients[2,1]))
P<-apply(data,1,function(x) summary(bayesglm(as.numeric(as.factor(anno$Gender))~as.numeric(x)))$coefficients[2,4])
Risk<-data.frame(OR,P)
Risk

setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/RnaseqBacterial/manuscript/result/strain")
data<-read.table("strain.tab.matrix.freq.txt",head=T,row.names=1,sep="\t",check.names=F)

per<-unlist(apply(data,1,function(x) sum(x>0)/length(x)))
ave<-round(unlist(apply(data,1,function(x) mean(x))),2)
sd<-round(unlist(apply(data,1,function(x) sd(x))),2)
iqr<-round(unlist(apply(data,1,function(x) IQR(x))),2)
min<-round(unlist(apply(data,1,function(x) min(x))),2)
max<-round(unlist(apply(data,1,function(x) max(x))),2)
out<-data.frame(per,ave,sd,iqr,min,max)
out<-out[order(out$ave,decreasing = T),]
out$per<-percent(out$per)
write.csv(out,file="strain.table.csv",quote=F,row.names = T)

sam<-read.csv("https://raw.githubusercontent.com/Shicheng-Guo/RnaseqBacterial/master/extdata/saminfo.csv")
colnames(data)<-unlist(lapply(strsplit(colnames(data),"_"), function(x) x[1]))
data<-data[,match(sam$SID,colnames(data))]
head(data)
data<-RowNARemove(data)
OR<-exp(apply(data,1,function(x) summary(bayesglm(as.numeric(as.factor(sam$ST))~as.numeric(x)))$coefficients[2,1]))
P<-apply(data,1,function(x) summary(bayesglm(as.numeric(as.factor(sam$ST))~as.numeric(x)))$coefficients[2,4])
Risk<-data.frame(OR,P)
Risk<-Risk[order(Risk$P,decreasing = F),]
out<-subset(Risk,P<0.05)
write.csv(out,file="strain.ST.sig.table.csv",quote=F,row.names = T)

OR<-exp(apply(data,1,function(x) summary(bayesglm(as.numeric(as.factor(sam$Severity))~as.numeric(x)))$coefficients[2,1]))
P<-apply(data,1,function(x) summary(bayesglm(as.numeric(as.factor(sam$Severity))~as.numeric(x)))$coefficients[2,4])
Risk<-data.frame(OR,P)
Risk<-Risk[order(Risk$P,decreasing = F),]
out<-subset(Risk,P<0.05)
write.csv(out,file="strain.severity.table.csv",quote=F,row.names = T)

