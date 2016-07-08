CharToNumberFactor<- function(Variable){ 
  Variable<-as.numeric(Variable)
  Variable[is.na(Variable)]<-0
  Variable
}
