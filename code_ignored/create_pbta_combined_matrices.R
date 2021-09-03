# Author: Komal S. Rathi
# Function: Create PBTA expected count and tpm matrix with unique gene symbols
# This needs to be run with every updated version of OpenPBTA (current v18)

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(dplyr))

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))
source(file.path(code_dir, 'code_ignored', 'collapse_matrix.R'))
pbta.dir <- file.path(ref_dir, 'pbta')

# adapt file to remove unwanted samples
pbta_adapt <- read.delim(file.path(pbta.dir, 'pbta-histologies-base-adapt.tsv'))
pbta_adapt <- pbta_adapt %>%
  filter(kf_visibility == "True")

# count matrix for polyA + stranded data
pbta.polya.counts <- readRDS(file.path(pbta.dir, 'pbta-gene-counts-rsem-expected_count.polya.rds'))
pbta.st.counts <- readRDS(file.path(pbta.dir, 'pbta-gene-counts-rsem-expected_count.stranded.rds'))

# collapse and merge
pbta.polya.counts <- collapse.rnaseq(pbta.polya.counts)
pbta.st.counts <- collapse.rnaseq(pbta.st.counts)
pbta.full.counts <- pbta.polya.counts %>%
  rownames_to_column('sym') %>%
  inner_join(pbta.st.counts %>%
               rownames_to_column('sym'), by = "sym") %>%
  column_to_rownames('sym')

# remove unwanted samples
pbta.full.counts <- pbta.full.counts[,colnames(pbta.full.counts) %in% pbta_adapt$Kids_First_Biospecimen_ID]

# save collapsed count matrix
saveRDS(pbta.full.counts, file = file.path(pbta.dir, 'pbta-gene-expression-rsem-counts-collapsed.polya.stranded.rds'))

# tpm matrix for polyA + stranded data
pbta.tpm.stranded <- readRDS(file.path(pbta.dir, 'pbta-gene-expression-rsem-tpm.stranded.rds'))
pbta.tpm.polya <- readRDS(file.path(pbta.dir, 'pbta-gene-expression-rsem-tpm.polya.rds'))

# collapse and merge
pbta.tpm.stranded <- collapse.rnaseq(pbta.tpm.stranded)
pbta.tpm.polya <- collapse.rnaseq(pbta.tpm.polya)

# remove unwanted samples
pbta.tpm.stranded <- pbta.tpm.stranded[,colnames(pbta.tpm.stranded) %in% pbta_adapt$Kids_First_Biospecimen_ID]
pbta.tpm.polya <- pbta.tpm.polya[,colnames(pbta.tpm.polya) %in% pbta_adapt$Kids_First_Biospecimen_ID]

# merge
pbta.full.tpm <- pbta.tpm.polya %>%
  rownames_to_column('sym') %>%
  inner_join(pbta.tpm.stranded %>%
               rownames_to_column('sym'), by = "sym") %>%
  column_to_rownames('sym')

# save
saveRDS(pbta.tpm.stranded, file = file.path(pbta.dir, 'pbta-gene-expression-rsem-tpm-collapsed.stranded.rds'))
saveRDS(pbta.tpm.polya, file = file.path(pbta.dir, 'pbta-gene-expression-rsem-tpm-collapsed.polya.rds'))
saveRDS(pbta.full.tpm, file = file.path(pbta.dir, 'pbta-gene-expression-rsem-tpm-collapsed.polya.stranded.rds'))

