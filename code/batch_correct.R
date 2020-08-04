# Author: Komal S. Rathi
# Function: Batch correct using ComBat

library(sva)

batch.correct <- function(mat, clin) {
  corrected.mat <- suppressWarnings(ComBat(dat = log2(mat + 1), batch = clin$batch))
  corrected.mat <- 2^(corrected.mat)
  return(corrected.mat)
}
