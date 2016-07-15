barplotReAxis<-function(mp){
mp<-c(0,mp)
np<-c()
end<-length(mp)-1
for(i in 1:end){
  np.tmp<-(mp[i]+mp[i+1])/2
  np<-c(np,np.tmp)
}
return(round(np))
}
