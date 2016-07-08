calcOddsRatio <- function(mymatrix,alpha=0.05,referencerow=2){
  numrow <- nrow(mymatrix)
  myrownames <- rownames(mymatrix)
  for (ii in 1:numrow){
    rowname <- myrownames[ii]
    DiseaseUnexposed <- mymatrix[referencerow,1]
    ControlUnexposed <- mymatrix[referencerow,2]
    if (ii != referencerow){ 
      DiseaseExposed <- mymatrix[ii,1]
      ControlExposed <- mymatrix[ii,2] 
      totExposed <- DiseaseExposed + ControlExposed
      totUnexposed <- DiseaseUnexposed + ControlUnexposed
      probDiseaseGivenExposed <- DiseaseExposed/totExposed
      probDiseaseGivenUnexposed <- DiseaseUnexposed/totUnexposed
      probControlGivenExposed <- ControlExposed/totExposed
      probControlGivenUnexposed <- ControlUnexposed/totUnexposed
      # calculate the odds ratio
      oddsRatio <- (probDiseaseGivenExposed*probControlGivenUnexposed)/(probControlGivenExposed*probDiseaseGivenUnexposed)
      # calculate a confidence interval
      confidenceLevel <- (1 - alpha)*100
      sigma <- sqrt((1/DiseaseExposed)+(1/ControlExposed)+(1/DiseaseUnexposed)+(1/ControlUnexposed))
      # sigma is the standard error of our estimate of the log of the odds ratio
      z <- qnorm(1-(alpha/2))
      lowervalue <- oddsRatio * exp(-z * sigma)
      uppervalue <- oddsRatio * exp( z * sigma)
      OR=paste(round(oddsRatio,2)," (",round(lowervalue,2),",",round(uppervalue,2),")",sep="")
    }
  }
  OR
}

