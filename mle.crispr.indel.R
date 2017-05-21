########################################################################################
###   Title: Estimate ture Indel freqency of CRSCPR in Case adjusted by Background Noise
###   Author: Shicheng Guo, Ph.D. Email: Shicheng.Guo@hotmail.com 
###   updata time: 9/1/2015
########################################################################################

mle<-function(n,R,q){
  #n=200   # reads number with indel  
  #R=1000  # total reads number 
  #q=0.10  # observed indel ratio in control 
  p=n/R   # observed indel ratio in case 
  n/R
  q11<-c()
  range<-seq(0,R,1)
  for(p0 in range){
    q1<-dbinom(n-p0,R-p0,q)
    q11=c(q11,q1)
  }
  plot(q11)
  mlefreq=range[which.max(q11)]/R
  z = 1.96
  ub = (R * mlefreq + z**2 / 2 + z * sqrt(R *mlefreq * (1 - mlefreq) + z**4 / 4)) / (R + z**2)
  lb = (R * mlefreq + z**2 / 2 - z * sqrt(R *mlefreq * (1 - mlefreq) + z**4 / 4)) / (R + z**2)
  return(c(mle=mlefreq,ub=ub,lb=lb))
}

# Example Usaage:
result<-mle(n=800,R=8000,q=0.05)
