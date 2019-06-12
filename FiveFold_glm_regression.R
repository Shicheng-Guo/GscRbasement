
input<-data
cv.error <- NULL
k <- 5
rlt1<-c()
rlt2<-c()
V<-c()
for(i in 1:k){
  index <- sample(1:nrow(input),round(0.9*nrow(input)))
  train.cv <- input[index,]
  test.cv <- input[-index,]
  fit <- (bayesglm(Phen ~ MethylHBV+Sex+Age,family=binomial,data=train.cv,na.action=na.omit))
  pscores <- predict(fit,test.cv)
  V=rbind(V,data.frame(test.cv$phen,pscores))
}
V
plotROC(data=V,cOutcome=1,predrisk=V$pscores)
