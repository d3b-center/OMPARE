# Author: Komal S. Rathi
# Function: Collapse matrix of gene ids to unique gene symbols
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(dplyr))

# function to collapse matrix
collapse.rnaseq <- function(expr.mat){
  # separate gene_id and gene_symbol
  expr.mat <- expr.mat %>% 
    mutate(gene_id = str_replace(gene_id, "_PAR_Y_", "_"))  %>%
    separate(gene_id, c("gene_id", "gene_symbol"), sep = "\\_", extra = "merge") %>%
    unique()
  
  # uniquify gene_symbol
  expr.collapsed <- expr.mat %>% 
    mutate(means = rowMeans(dplyr::select(.,-gene_id, -gene_symbol))) %>% # take rowMeans
    arrange(desc(means)) %>% # arrange decreasing by means
    distinct(gene_symbol, .keep_all = TRUE) %>% # keep the ones with greatest mean value. If ties occur, keep the first occurencce
    dplyr::select(-c(means, gene_id)) %>%
    unique() %>%
    remove_rownames() %>%
    column_to_rownames('gene_symbol')
  
  return(expr.collapsed)
}
