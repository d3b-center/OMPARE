get_zscore <- function(x) {
  x <- log2(x+1)
  out <- (x-mean(x))/sd(x)
  return(out)
}
