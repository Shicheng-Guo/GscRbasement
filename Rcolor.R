redscale<-function(n=24){
  rgb((1:n) - 1, 0, 0, maxColorValue=n)
}

blueyellow<-function(n=24){
  color=colorRampPalette(c("blue", "yellow"))(n) 
  return(color)  
}

whitered<-function(n=24){
  color=colorRampPalette(c("white", "red"))(n) 
  return(color)  
}

whiteyellow<-function(n=24){
  color=colorRampPalette(c("white", "yellow"))(n) 
  return(color)  
}

whiteblue<-function(n=24){
  color=colorRampPalette(c("white", "blue"))(n) 
  return(color)  
}

whitebrown<-function(n=24){
  color=colorRampPalette(c("white", "brown"))(n) 
  return(color)  
}

whiteblack<-function(n=24){
  color=colorRampPalette(c("white", "black"))(n) 
  return(color)  
}

yellowblue<-function(n=24){
  color=colorRampPalette(c("yellow", "blue"))(n) 
  return(color)  
}

jetColors<-function(N=24){
  k <- ceiling(N/4)
  temp.red <- c(rep(0, 2 * k), 1:k, rep(k, k - 1), k:1)
  temp.green <- c(rep(0, k), 1:k, rep(k, k - 1), k:1, rep(0,k))
  temp.blue <- c(1:k, rep(k, k - 1), k:1, rep(0, 2 * k))
  temp.rgb <- cbind(temp.red, temp.green, temp.blue)
  delta <- 5 * k - 1 - N
  delta <- ceiling(delta/2)
  temp.rgb <- temp.rgb[delta:(delta + N - 1), ]/k
  color<-rgb(temp.rgb[, 1], temp.rgb[, 2], temp.rgb[, 3])
  return(color)
}
