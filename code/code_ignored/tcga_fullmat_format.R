###########################
# Format TCGA data for GSVA
###########################

library(tidyverse)
library(RDiseaseXpress)

# TCGA TPM data
dat <- readRDS('data/Reference/TCGA/tumor_datasets_TPM.RDS')
tcga.meta <- RDiseaseXpress::getSamples(myStudy = 'TCGA')
tcga.meta <- tcga.meta[grep('Primary', tcga.meta$definition),]
tcga.tpm <- dat[,colnames(dat) %in% tcga.meta$sample_id]
tcga.meta <- tcga.meta %>%
  mutate(library_name = "polyA") %>%
  filter(tcga.meta$sample_id %in% colnames(tcga.tpm)) %>%
  dplyr::select(sample_id, sample_barcode, ethnicity, race, vital_status, gender, disease, disease_name, disease_subtype, overall_survival_time_in_days, library_name) %>%
  column_to_rownames('sample_id')
tcga.tpm <-  tcga.tpm[,rownames(tcga.meta)]
saveRDS(tcga.tpm, file = 'data/Reference/TCGA/TCGA_matrix_TPM.RDS')
saveRDS(tcga.meta, file = 'data/Reference/TCGA/TCGA_meta.RDS')

###########################
# Format TCGA data for TIS
###########################

# TCGA count data
tcga.counts <- readRDS('~/Projects/toil-rnaseq-20k/data/results/expression/tumor_datasets_counts.RDS')
tcga.counts <-  tcga.counts[,rownames(tcga.meta)]
saveRDS(tcga.counts, file = 'data/Reference/TCGA/TCGA_matrix_counts.RDS')
