# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", 'quiet.R'))

# xcell deconvolution
immune_profile <- function(fullmat) {
  # immune scores
  xcell_scores <- as.data.frame(quiet(xCellAnalysis(expr = fullmat)))
  return(xcell_scores)
}
