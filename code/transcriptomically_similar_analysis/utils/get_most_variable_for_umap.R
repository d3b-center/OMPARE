# get 10000 most variable genes for umap plotting & nearest neighbor analysis
get_most_variable_for_umap <-  function(expr_corrected){
  # filter out low expression genes
  maxVals <- apply(expr_corrected, FUN = max, MARGIN = 1)
  expr_filtered <- expr_corrected[maxVals > 20,]
  
  # for dimensionality reduction visualization (dim_reduction_plot.R)
  rv <- matrixStats::rowVars(as.matrix(expr_filtered))
  select <- order(rv, decreasing=TRUE)[seq_len(10000)]
  expr_most_var <- expr_filtered[select,]
  
  # return output
  return(expr_most_var)
}