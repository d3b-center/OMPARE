# function for upper quartile (UQ) normalization function
UQ <- function(X){
  uq<-function(y){
    quantile(y, 0.75)
  }
  X<-X+0.1
  upperQ<-apply(X,2,uq)
  f<-upperQ/mean(upperQ) # calculate normalization factor
  res<-scale(X,center=FALSE,scale=f) 
  return(res)
}
