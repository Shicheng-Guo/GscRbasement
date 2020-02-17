#install.packages("quantmod")
#install.packages("TTR")
library(quantmod)
library(TTR)
myATR        <- function(x) ATR(HLC(x))[,'atr']
mySMI        <- function(x) SMI(HLC(x))[, "SMI"]
myADX        <- function(x) ADX(HLC(x))[,'ADX']
myAroon      <- function(x) aroon(cbind(Hi(x),Lo(x)))$oscillator
myBB         <- function(x) BBands(HLC(x))
myChaikinVol <- function(x) Delt(chaikinVolatility(cbind(Hi(x),Lo(x))))[, 1]
myEMA        <- function(x,n) EMA(HLC(x)[,3],n)
myEMV        <- function(x) EMV(cbind(Hi(x),Lo(x)),Vo(x))[,2]
myMACD       <- function(x) MACD(Cl(x),percent=F,nFast=12, nSlow=26, nSig=9)
myMFI        <- function(x) MFI(HLC(x), Vo(x))
mySAR        <- function(x) SAR(cbind(Hi(x),Cl(x))) [,1]
myVolat      <- function(x) volatility(OHLC(x),calc="garman")[,1]
myCCI        <- function(x) CCI(x,n = 20, maType, c = 0.015)
myRSI        <- function(x) RSI(x[,4], n=2)

CrossOver<-function(x,y,n=2){
  return(tail(x,1)[1,1]>=tail(y,1)[1,1] && tail(x,n)[1,1]<=tail(y,n)[1,1])
}

beginning <- as.Date(Sys.Date()-450)
date<-Sys.Date()

#setwd("/home/sguo/Dropbox/stock/mystock")
setwd("N:\\R")
input<-read.table("bak.db",sep="\t")

#setwd("C:\\Users\\shg047\\Dropbox\\stock\\mystock")
#input<-read.table("C:\\Users\\shg047\\Dropbox\\stock\\bak.db",sep="\t")
rlt<-c()
for(i in 1:nrow(input)){
  symbol=as.character(input[i,1])
  try(Object<-getSymbols(Symbols=as.character(symbol),src="google",env=NULL,from=beginning,to=date,auto.assign=FALSE),silent = T)
  #if(nrow(na.omit(Object))<300) next
  #if(any(is.na(tail(Object,1)[1,1:3]))) next
  if((tail(Object,1)[1,5])<50000) next
  if((nrow(Object)<60)) next
  #if((tail(Object,1)[1,1])>500) next
  #if((tail(Object,1)[1,1])<5) next
  v.now<-mean(na.omit(head(tail(Object,15),15)[,5]))
  v.before<-mean(na.omit(head(tail(Object,60),45)[,5]))
  p.now<-mean(na.omit(head(tail(Object,15),15)[,1]))
  p.before<-mean(na.omit(head(tail(Object,60),45)[,1]))
  
  if(v.now/v.before>10 && p.now>1.10*p.before){
    print(symbol)
    #file=paste(symbol,"png",sep=".")
    #png(file,width = 8, height = 8, units = 'in', res = 800)
    #chartSeries(Object,subset='last 6 months',main=symbol,TA=c(addVo(),addMACD(),addBBands(),addCCI(),addRSI(n = 12, wilder = TRUE),addEMA(n=60, col="blue"),addEMA(n=120, col="yellow")))
    #dev.off()
  }
  #myEMAD200<-myEMA(na.omit(Object),n=200)
  #myEMAD5<-myEMA(na.omit(Object),n=5)
  #myEMAD10<-myEMA(na.omit(Object),n=10)
  #myEMAD20<-myEMA(na.omit(Object),n=20)
  #myEMAD30<-myEMA(na.omit(Object),n=30)
  #myEMAD60<-myEMA(na.omit(Object),n=60)
  #myEMAD120<-myEMA(na.omit(Object),n=120)
  #myEMAD200<-myEMA(na.omit(Object),n=200)
  #if(tail(Object,1)[1,4]<tail(myEMAD200,1)[1,1]) next
  #if(tail(Object,1)[1,4]<tail(myEMAD120,1)[1,1]) next
  #if(tail(Object,1)[1,4]<tail(myEMAD120,1)[1,1]) next
  #if(tail(Object,1)[1,4]<tail(myEMAD60,1)[1,1]) next
  #if(! all(tail(myEMAD5,10)[,1]>tail(myEMAD10,10)[,1])) next
  #if(! all(tail(myEMAD10,10)[,1]>tail(myEMAD30,10)[,1])) next
  #if(! all(tail(myEMAD30,10)[,1]>tail(myEMAD60,10)[,1])) next
  #BBobject<-myBB(na.omit(Object))         
  #if(mean(tail(BBobject,2)[,3]-tail(BBobject,2)[,1])/mean(tail(BBobject,30)[,3]-tail(BBobject,30)[,1])<1.05) next
  #if(! tail(myEMAD5,1)[1,1]>=tail(myEMAD10,1)[1,1] && tail(myEMAD5,5)[1,1]<=tail(myEMAD10,5)[1,1]) next
  #if(mean(tail(BBobject,2)[,3]-tail(BBobject,2)[,1])/mean(tail(BBobject,30)[,3]-tail(BBobject,30)[,1])>1.20) next
  #if(! tail(myEMAD10,1)[1,1]>=tail(myEMAD20,1)[1,1] && tail(myEMAD10,5)[1,1]<=tail(myEMAD20,5)[1,1]) next
  #if( ! (CrossOver(myEMAD5,myEMAD20) && CrossOver(myEMAD10,myEMAD20) && CrossOver(myEMAD20,myEMAD60))) next
  #if(CrossOver(myEMAD5,myEMAD20)) rlt<-rbind(rlt,c(symbol,"5cross20"))
  #if(CrossOver(myEMAD5,myEMAD30)) rlt<-rbind(rlt,c(symbol,"5cros30"))
  #if(CrossOver(myEMAD5,myEMAD60)){
  #  rlt<-rbind(rlt,c(symbol,"5cross60"))
  #  print(c(i,symbol))
  #  file=paste(symbol,"png",sep=".")
  #  png(file,width = 8, height = 8, units = 'in', res = 800)
  #  chartSeries(Object,subset='last 12 months',main=symbol,TA=c(addVo(),addMACD(),addBBands(),addCCI(),addRSI(n = 12, wilder = TRUE),addEMA(n=60, col="blue"),addEMA(n=120, col="yellow")))
  #  dev.off()
  #}
  
  #if(CrossOver(myEMAD5,myEMAD120)) rlt<-rbind(rlt,c(symbol,"5cross120"))
  #if(CrossOver(myEMAD5,myEMAD200)) rlt<-rbind(rlt,c(symbol,"5cross200"))
  #if(CrossOver(myEMAD5,myEMAD60)) print(c(symbol,"5cross60"))
  
  # daily to weekly or monthly
  # Object1<-to.weekly(Object)
  # Object1<-to.monthly(Object)
  #MACDobject<-myMACD(Cl(Object1))    
  #cond1<- as.numeric(MACDobject[nrow(MACDobject),1]) >= -1 && as.numeric(MACDobject[nrow(MACDobject)-2,1]) < -2
  # MSIDobject<-mySMI(na.omit(Object))
  #BBobject<-myBB(na.omit(Object))   
  #tmp<-tail(BBobject,n=8)[,3]-tail(BBobject,n=8)[,1]
  #if(sd(tmp[1:8])/sd(tmp[1:7])<2) next
  #if(tmp[8]/mean(tmp[1:7],na.rm=T)>1.2) next
  
  #print(c(i,symbol))
  #CCIobject<-CCI(na.omit(HLC(Object1)))
  #myRSIValue<-myRSI(na.omit(Object1))
  #tail(CCIobject)
  #tail(myRSIValue)
  #if(Object[nrow(Object),5]>100000 && cond7>=1.1 && cond7<=1.5){
  #if(Object[nrow(Object),5]>100000 && tail(CCIobject,n=1)[1]< -95 && tail(myRSIValue,n=1)[1]<45){
} 
#file=paste("A5up60",date,"txt",sep=".")
#filewrite.table(rlt,file,quote=F,sep="\t",row.names = F,col.names = F)


#beginning <-"2017-07-01"
#date<-"2017-08-30"
#beginning <-"2016-10-11"
#date<-"2016-12-11"
#symbol="QBAK"
#try(Object<-getSymbols(Symbols=as.character(symbol),src="google",env=NULL,from=beginning,to=date,auto.assign=FALSE),silent = T)
#tail(Object)
#BBobject<-myBB(na.omit(Object))
#head(BBobject)
#tail(BBobject)
#file=paste(symbol,"png",sep=".")
#png(file,width = 8, height = 8, units = 'in', res = 800)
#chartSeries(Object,subset='last 6 months',main=symbol,TA=c(addVo(),addMACD(),addBBands(),addCCI()))
#dev.off()



