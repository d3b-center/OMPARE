# Author: Komal S. Rathi
# Function: 
# 1. Create PBTA expected count and tpm matrix with unique gene symbols

library(tidyverse)
library(dplyr)
source('code/code_ignored/collapse_matrix.R')

# count matrix for polyA + stranded data
pbta.polya.counts <- readRDS('data/Reference/PBTA/pbta-gene-counts-rsem-expected_count.polya.rds')
pbta.st.counts <- readRDS('data/Reference/PBTA/pbta-gene-counts-rsem-expected_count.stranded.rds')

# collapse and merge
pbta.polya.counts <- collapse.rnaseq(pbta.polya.counts)
pbta.st.counts <- collapse.rnaseq(pbta.st.counts)
pbta.full.counts <- pbta.polya.counts %>%
  rownames_to_column('sym') %>%
  inner_join(pbta.st.counts %>%
               rownames_to_column('sym'), by = "sym") %>%
  column_to_rownames('sym')

# save collapsed count matrix
saveRDS(pbta.full.counts, file = 'data/Reference/PBTA/pbta-gene-expression-rsem-counts-collapsed.polya.stranded.rds')

# tpm matrix for polyA + stranded data
pbta.tpm.stranded <- readRDS('data/Reference/PBTA/pbta-gene-expression-rsem-tpm.stranded.rds')
pbta.tpm.polya <- readRDS('data/Reference/PBTA/pbta-gene-expression-rsem-tpm.polya.rds')

# collapse and merge
pbta.tpm.stranded <- collapse.rnaseq(pbta.tpm.stranded)
pbta.tpm.polya <- collapse.rnaseq(pbta.tpm.polya)

# merge
pbta.full.tpm <- pbta.tpm.polya %>%
  rownames_to_column('sym') %>%
  inner_join(pbta.tpm.stranded %>%
               rownames_to_column('sym'), by = "sym") %>%
  column_to_rownames('sym')

# save
saveRDS(pbta.tpm.stranded, file = 'data/Reference/PBTA/pbta-gene-expression-rsem-tpm-collapsed.stranded.rds')
saveRDS(pbta.tpm.polya, file = 'data/Reference/PBTA/pbta-gene-expression-rsem-tpm-collapsed.polya.rds')
saveRDS(pbta.full.tpm, file = 'data/Reference/PBTA/pbta-gene-expression-rsem-tpm-collapsed.polya.stranded.rds')

