# Author: Komal S. Rathi
# Function: Batch correct using ComBat

library(sva)

batch.correct <- function(mat, clin) {
  if(identical(rownames(clin), colnames(mat))){
    print("Matching dimensions")
  } else {
    print("Check inputs")
    break
  }
  
  corrected.mat <- ComBat(dat = log2(mat + 1), batch = clin$batch)
  corrected.mat <- 2^(corrected.mat)
  return(corrected.mat)
}
