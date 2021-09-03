# extract nearest neighbor info
extract_umap_nearest_neighbor_table <- function(umap_out, expr_most_var, patient_of_interest){
  corr <- as.data.frame(umap_out$nn$correlation$idx) # nn
  dist <- as.data.frame(umap_out$nn$correlation$dist) # distances
  corr <- t(apply(corr, MARGIN = 1, FUN = function(x) colnames(expr_most_var)[x]))
  nn_table <- data.frame(nearest_neighbor = as.character(corr[grep(patient_of_interest, rownames(corr)),]), 
                         distance = as.numeric(dist[grep(patient_of_interest, rownames(dist)),]))
  nn_table$distance <- round(nn_table$distance, digits = 3)
  
  # nearest neighbor table
  return(nn_table)
}

