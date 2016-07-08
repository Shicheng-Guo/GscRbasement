frr<-function(rr1,rr2,p,q,lamdi,frr){
#familial reltaive risk
r1<-10   #relative risk(estimated by odds ratios) for hetrozygotes relative to common homozygotes
r2<-10   #relative risk(estimated by odds ratios) for rare homozygotes to common homozygotes
p<-0.3   #risk allelefrequency
q<-1-p
lamdi0<-8.48  #overall familial relatve risk 
frr<-(p*(p*r2+q*r1)^2+q*(p*r1+q*r2)^2)/(p^2*r2+2*p*q*r1+q^2)^2
proportion<-log(lamdi)/log(lamdi0)
frr
proportion
}
