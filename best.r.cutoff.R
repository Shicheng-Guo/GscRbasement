# estimate the best cor cutoff threshold for findCorrelation() for EPIC dataset (Atrial Fibrillation)

input<-t(myNorm[na.omit(match(subset(map,map[1]=="chr22")[,4],rownames(myNorm))),])
R<-c()
P<-c()
for(i in 1:ncol(input)){
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
pdf("P-R.pdf")
plot(-log(P,10),R)
dev.off()
