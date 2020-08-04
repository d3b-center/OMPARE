quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

getZ <- function(x) {
  x <- log2(x+1)
  out <- (x-mean(x))/sd(x)
  return(out)
}