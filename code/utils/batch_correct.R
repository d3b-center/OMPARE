# Author: Komal S. Rathi
# Function: Batch correct using ComBat

suppressPackageStartupMessages({
  library(sva)
})

batch_correct <- function(mat, clin) {
  corrected_mat <- suppressWarnings(ComBat(dat = log2(mat + 1), batch = clin$batch))
  corrected_mat <- 2^(corrected_mat)
  return(corrected_mat)
}
