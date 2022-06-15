# Function to calculate v-scores for gene and cluster 

suppressPackageStartupMessages({
  library(tidyverse)
})

####### Define function to calculate v-score
# adapted from https://github.com/d3b-center/
# d3b-pnoc003-HGG-DMG-omics/blob/master/analyses/immune-deconv/util/calc_vtest.R)
# function to compute v.test
compute.v.test <- function(x, clustering_col){
  
  # get unique clusters i.e. 1-n clusters
  clusters <- unique(x[[clustering_col]])
  
  out <- data.frame()
  # iterate over each cluster
  for(i in 1:length(clusters)){
    
    # subset x to cluster i
    y = x %>% 
      filter(get(clustering_col) == clusters[i])
    
    # mean gene expression per cluster
    mean.gene.cluster = unique(y$cluster_gene_mean_score)
    
    # mean gene expression (i.e. global mean)
    mean.gene = unique(y$gene_mean_score)
    
    # calculate numerator
    num = mean.gene.cluster - mean.gene
    
    # total sample size
    n = nrow(x)
    
    # cluster sample size
    ng = nrow(y)
    
    # variance of gene expression (i.e. global variance)
    var.gene = unique(x$gene_variance) 
    
    # calculate denominator
    denom = (n-ng/n-1)*(var.gene/ng)
    denom = sqrt(denom)
    
    # calculate vscore
    v = num/denom
    out[i,'cluster'] <- clusters[i]
    out[i,'v_score'] <- v
  }
  return(out)
}
