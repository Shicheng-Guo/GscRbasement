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

beginning <- as.Date(Sys.Date()-4000)
date<-Sys.Date()
symbol=as.character("VXX")
try(object<-getSymbols(Symbols=as.character(symbol),src="yahoo",env=NULL,from=beginning,to=date,auto.assign=FALSE),silent = T)

tail(object)

par(mfrow=c(3,1))

for(day in c(90,60,30)){
c90<-c()
for(i in 1:(length(Object[,1])-day)){
  c<-as.numeric(Object[i+day,4])/as.numeric(Object[i,4])
  c90<-c(c90,c)
}
barplot(c90-1)
}

vector<-c90-1
vector<-vector[vector<0]
quantile(vector)





















hist(c90-1)
median(c90-1)
mean(c90-1)

vix<-Object[,4]-Object[,1]>0
vix30<-c()

for(i in 1:(length(vix)-30)){
  vi<-sum(vix[i:(i+30)])
  vix30<-c(vix30,vi)
}
plot(Object[,4],pch=16,cex=0.5)

quantile(Object[,4],seq(0,1,0.1))

pop<-which(Object[,4]>30)
Object[pop,]
tail(pop)
POP<-pop[1:(length(pop)-3)]
tail(POP)

V<-c()
for(i in seq(10,62,3)){
  S=as.numeric(Object[POP,4])
  E=as.numeric(Object[POP+i,4])
  V=cbind(V,(S-E)/S)
}
colnames(V)=seq(10,62,3)
head(V)
boxplot(V)










