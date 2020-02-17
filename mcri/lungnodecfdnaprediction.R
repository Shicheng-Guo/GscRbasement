setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/lungnodes/extdata")

data<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/lungnodes/master/extdata/TCGA.3108.Marker.DJG.beta.txt",head=T)
colnames(data)[2]<-"phen"
iid<-data[,1]
phen<-data[,2]
table(phen)
clinical<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/lungnodes/master/extdata/clinical.txt",head=T,row.names = 1)
input<-data.matrix(data[,2:ncol(data)])
rownames(input)<-data[,1]
input[1:5,1:5]
input[input[,1]<3,1]<-0
input[input[,1]==3,1]<-1
input[1:5,1:5]
mean(input[,1])
input<-input[,-which(apply(input,2,function(x) mean(x,na.rm=T))>0.61)]
input[,1]
match(iid,rownames(clinical))

rlt<-c()
for(i in 2:ncol(input)){
  fit<-summary(glm(as.factor(input[,1])~input[,i],family=binomial))$coefficients[2,]
  rlt<-rbind(rlt,fit)
  print(i)
}
rlt<-data.frame(rlt)
rownames(rlt)<-colnames(input)[2:ncol(input)]
rlt<-subset(rlt,Estimate>0)
rlt<-rlt[order(rlt[,4]),]
head(rlt)

temp<-input[,c(1,match(rownames(rlt[1:150,]),colnames(input)))]
write.csv(temp,file="lungnode.scores.csv")

temp[1:5,1:5]
rownames(temp)<-temp[,1]
temp[1:5,1:5]
table(temp[,1])
head(temp)

heatmap((temp[,2:ncol(temp)]))
xx<-apply(temp[,2:ncol(temp)],1,sum)
xx<-xx[order(names(xx))]
plot(xx,col=as.numeric(names(xx))+1)

address<-lapply(strsplit(colnames(input)[2:ncol(input)],"_"),function(x) paste(x[1],":",x[2],"-",as.numeric(x[2])+1,sep=""))
write.table(t(data.frame(address)),file="biomarker3108.txt",sep="\t",quote=F,col.names = F,row.names = F)

#### find suspect samples
library("tidyverse")
library("caret")
library("PredictABEL")
library("MASS")

temp<-input[,c(1,match(rownames(rlt[1:250,]),colnames(input)))]
temp[1:5,1:5]
xx<-apply(temp[,2:ncol(temp)],1,sum)
matrix<-data.frame(score=xx,phen=temp[,1])
plot(matrix$score,col=matrix$phen+1,type="n")
text(rownames(matrix),x=1:nrow(matrix),y=matrix[,1],cex=0.75,col=matrix$phen+1)

subset(matrix,score>5.5 & phen==0)
subset(matrix,score<4.5 & phen==1)

glm.fit<-glm(as.factor(phen)~.,data=data.frame(temp),family=binomial)
step.model <- glm.fit %>% stepAIC(trace = FALSE)
coef(step.model)
pred1 <- predRisk(step.model)
par(mfrow=c(2,2),cex.lab=1.5,cex.axis=1.5)
plotROC(data=temp,cOutcome=1,predrisk=cbind(pred1))

matrix<-data.frame(phen=temp[,1],predict=pred1)
head(matrix)
head(xx)
plot(matrix[,1]~matrix[,2])
#############

index=data.frame(iid=rownames(temp),phen=temp[,1],score=apply(temp[,2:ncol(temp)],1,sum))
write.table(index,file="lungnodes.index.txt",sep="\t",col.names = T,row.names = F,quote=F)             
getwd()             

input<-input[,c(1,match(rownames(rlt[1:30,]),colnames(input)))]
dim(input)
input<-input[,-match("chr2_74153317",colnames(input))]
rownames(input)<-input[,1]
heatmap(t(input[,2:ncol(input)]))

plot(density(input[which(input[,1]==0),2]),col="blue",xlim=c(0,0.15))
lines(density(input[which(input[,1]==1),2]),col="red",xlim=c(0,0.15))

library("randomForest")
library("neuralnet")
library("arm")
library("PredictABEL")

set.seed(49)
cv.error <- NULL
k <- 2
rlt1<-c()
rlt2<-c()
for(i in 1:k){
  index <- sample(1:nrow(input),round(0.9*nrow(input)))
  train.cv <- input[index,]
  test.cv <- input[-index,]
  train.cv[1:3,1:3]
  test.cv[1:3,1:3]
  
  P=apply(train.cv[,2:ncol(train.cv)],2,function(x) summary(glm(as.factor(train.cv[,1])~x,family=binomial))$coefficients[2,4])
  
  train.cv<-train.cv[,c(1,na.omit(match(names(P[head(order(P),n=5200)]),colnames(train.cv))))]
  test.cv<-test.cv[,c(1,na.omit(match(names(P[head(order(P),n=5200)]),colnames(test.cv))))]

  train.cv[1:5,]
  
  RF <- randomForest(as.factor(phen) ~ ., data=data.frame(na.omit(train.cv)), importance=TRUE,proximity=T)
  imp<-RF$importance
  head(imp)
  imp<-imp[order(imp[,4],decreasing = T),]
  
  if(nrow(imp)>=30){
  topvar<-match(rownames(imp)[1:30],colnames(input))
  }else{
  topvar<-match(rownames(imp),colnames(input))
  }
  train.cv <- input[index,c(1,topvar)]
  test.cv <- input[-index,c(1,topvar)]
  
  RF <- randomForest(as.factor(phen) ~ ., data=data.frame(na.omit(train.cv)), importance=TRUE,proximity=T)
  
  n <- colnames(train.cv)
  f <- as.formula(paste("phen ~", paste(n[!n %in% "phen"], collapse = " + ")))
  
  nn <- neuralnet(f,data=train.cv,hidden=c(5,3),act.fct = "logistic",linear.output = T)
  pr.nn <- neuralnet::compute(nn,test.cv)
  trainRlt<-data.frame(phen=train.cv[,1],pred=unlist(nn$net.result[[1]][,1]))
  testRlt<-data.frame(phen=test.cv[,1],pred=unlist(pr.nn$net.result[,1]))
  rownames(trainRlt)=row.names(train.cv)
  rownames(testRlt)=row.names(test.cv)
  rlt1<-rbind(rlt1,trainRlt)  
  rlt2<-rbind(rlt2,testRlt)
  print(i)
}

data1<-na.omit(data.frame(rlt1))
data2<-na.omit(data.frame(rlt2))
model.glm1 <- bayesglm(phen~.,data=rlt1,family=binomial(),na.action=na.omit)
model.glm2 <- bayesglm(phen~.,data=rlt2,family=binomial(),na.action=na.omit)
pred1 <- predRisk(model.glm1)
pred2 <- predRisk(model.glm2)
par(mfrow=c(2,2),cex.lab=1.5,cex.axis=1.5)
plotROC(data=data1,cOutcome=1,predrisk=cbind(pred1))
plotROC(data=data2,cOutcome=1,predrisk=cbind(pred2))

### ROC 
pdf("mRNA.drugresponse.pdf")
par(mfrow=c(2,2),cex.lab=1.5,cex.axis=1.5)
plotROC(data=data1,cOutcome=1,predrisk=cbind(pred1))
plotROC(data=data2,cOutcome=1,predrisk=cbind(pred2))
dev.off()

### heatmap
library("gplots")
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/HeatMap.R")
input<-na.omit(input)
temp<-t(input[,2:ncol(input)])
colnames(temp)<-input[,1]
pdf("heatmap.lungnodes.pdf")
HeatMap(temp)
dev.off()

chr<-unlist(lapply(strsplit(colnames(input)[2:ncol(input)],"_"), function(x) x[1]))
pos<-unlist(lapply(strsplit(colnames(input)[2:ncol(input)],"_"), function(x) x[2]))
probe<-data.frame(CHR=chr,POS=pos)
head(probe)

heatmap(input[,2:ncol(input)])
