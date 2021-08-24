# xcell deconvolution

source(file.path(patient_level_analyses_utils, 'quiet.R'))

immune_profile <- function(fullmat) {
  # immune scores
  xcell_scores <- as.data.frame(quiet(xCellAnalysis(expr = fullmat)))
  return(xcell_scores)
}
