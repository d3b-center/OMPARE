# get umap output
get_umap_output <- function(expr_most_var){
  set.seed(100)
  umap_out <- uwot::umap(X = t(log2(expr_most_var+1)), n_neighbors = 21, n_components = 2, metric = "correlation", ret_nn = TRUE, n_sgd_threads = 123L)
  
  # add colnames/rownames to embeddings
  colnames(umap_out$embedding) <- c("UMAP1", "UMAP2")
  rownames(umap_out$embedding) <- colnames(expr_most_var)
  
  # add rownames to nearest neighbor
  rownames(umap_out$nn$correlation$idx) <- colnames(expr_most_var)
  rownames(umap_out$nn$correlation$dist) <- colnames(expr_most_var)
  
  # return output
  return(umap_out)
}