library("PredictABEL")

model.glm0 <- (bayesglm(phen ~ LINE.1,family=binomial,data=data,na.action=na.omit))
model.glm1 <- (bayesglm(phen ~ cg12071888,family=binomial,data=data,na.action=na.omit))
model.glm2 <- (bayesglm(phen ~ AHRR,family=binomial,data=data,na.action=na.omit))
model.glm3 <- (bayesglm(phen ~ ARRB2,family=binomial,data=data,na.action=na.omit))
model.glm4 <- (bayesglm(phen ~ FGF8,family=binomial,data=data,na.action=na.omit))
model.glm5 <- (bayesglm(phen ~ cg12071888+AHRR+ARRB2+FGF8,family=binomial,data=data,na.action=na.omit))

pred0 <- predRisk(model.glm0)
pred1 <- predRisk(model.glm1)
pred2 <- predRisk(model.glm2)
pred3 <- predRisk(model.glm3)
pred4 <- predRisk(model.glm4)
pred5 <- predRisk(model.glm5)

par(cex.lab=1.5,cex.axis=1.5,lwd=1.5)
plotROC(data=data,cOutcome=2,predrisk=cbind(pred0,pred1,pred2,pred3,pred4,pred5))
legend.text<-c("LINE-1","BARHL2","AHRR","ARRB2","FGF8","BAAF")
legend("bottomright",cex=1.5,legend=legend.text,lty=1:6,col=1:6,bty="n")

tiff("multi-AUC-splited.tiff")
pdf("multi-AUC-splited.pdf")
par(mfrow=c(2,3),cex.lab=1.5,cex.axis=1.5,lwd=1.5)
plotROC(data=data,cOutcome=2,predrisk=cbind(pred0),plottitle="LINE-1")
plotROC(data=data,cOutcome=2,predrisk=cbind(pred1),plottitle="BARHL2")
plotROC(data=data,cOutcome=2,predrisk=cbind(pred2),plottitle="AHRR")
plotROC(data=data,cOutcome=2,predrisk=cbind(pred3),plottitle="ARRB2")
plotROC(data=data,cOutcome=2,predrisk=cbind(pred4),plottitle="FGF8")
plotROC(data=data,cOutcome=2,predrisk=cbind(pred5),plottitle="BAAF")
dev.off()

predicted.value = predict(model.glm5,type=c("response"))
predicted.data  = data.frame(Type=data$phen, predicted.value)
logistic.rocobj  = roc(predicted.data$Type, predicted.data$predicted.value,smooth = FALSE)
logistic.rocdata = data.frame(Sens = logistic.rocobj$sensitivities, Spec = logistic.rocobj$specificities)
AUC= logistic.rocobj$auc[[1]]
logistic.rocdata[,3] = logistic.rocdata[,1] + logistic.rocdata[,2]
seq.max = which(logistic.rocdata[,3] == max(logistic.rocdata[,3]))
Sens= logistic.rocdata[seq.max[1],1]
Spec= logistic.rocdata[seq.max[1],2]
ACC=c(Sens,Spec, AUC)
ACC
