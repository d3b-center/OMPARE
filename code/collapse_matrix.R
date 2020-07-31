# Author: Komal S. Rathi
# Function: Collapse PBTA data to unique gene symbols
library(tidyverse)

# function to collapse matrix
collapse.rnaseq <- function(expr.mat){
  # separate gene_id and gene_symbol
  expr.mat <- expr.mat %>% 
    mutate(gene_id = str_replace(gene_id, "_PAR_Y_", "_"))  %>%
    separate(gene_id, c("gene_id", "gene_symbol"), sep = "\\_", extra = "merge") %>%
    unique()
  
  # uniquify gene_symbol
  expr.collapsed <- expr.mat %>% 
    mutate(means = rowMeans(select(.,-gene_id, -gene_symbol))) %>% # take rowMeans
    arrange(desc(means)) %>% # arrange decreasing by means
    distinct(gene_symbol, .keep_all = TRUE) %>% # keep the ones with greatest mean value. If ties occur, keep the first occurencce
    select(-c(means, gene_id)) %>%
    unique() %>%
    remove_rownames() %>%
    column_to_rownames('gene_symbol')
  
  return(expr.collapsed)
}

# PBTA
pbta.tpm.stranded <- readRDS('data/Reference/PBTA/pbta-gene-expression-rsem-tpm.stranded.rds')
pbta.tpm.stranded <- collapse.rnaseq(pbta.tpm.stranded)
saveRDS(pbta.tpm.stranded, file = 'data/Reference/PBTA/pbta-gene-expression-rsem-tpm-collapsed.stranded.rds')

# polyA
pbta.tpm.polya <- readRDS('data/Reference/PBTA/pbta-gene-expression-rsem-tpm.polya.rds')
pbta.tpm.polya <- collapse.rnaseq(pbta.tpm.polya)
saveRDS(pbta.tpm.polya, file = 'data/Reference/PBTA/pbta-gene-expression-rsem-tpm-collapsed.polya.rds')

