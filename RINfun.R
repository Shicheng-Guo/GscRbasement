# reverse-normal transform

RINfun=function(yorig){
  yranks=rank(yorig)
  tempp=(yranks-.5)/(length(yranks))
  return(qnorm(tempp))
}
