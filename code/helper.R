################################
# Miscellaneous helper functions
################################

# misc functions
getEns <- function(x) {
  out <- strsplit(x, split="_")[[1]][1]
  return(out)
}

remDotStuff <- function(x) {
  out <- strsplit(x, "\\.")[[1]][1]
}
