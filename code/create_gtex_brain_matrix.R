# Script to create GTEx Brain TPM matrix
library(tidyverse)

# GTEx Brain meta file
gtexClin <- readRDS('data/Reference/GTEx/GTEx_clinical.RDS')
gtexBrain <- gtexClin %>%
  filter(subtissue == "Brain") %>%
  mutate(tmp = sample_id)  %>%
  column_to_rownames('tmp')

# TPM matrix for all GTEx data
gtexData <- readRDS('data/Reference/GTEx/gtex_normals_TPM.RDS')
gtexData <- gtexData[,rownames(gtexBrain)]

# save matrices
saveRDS(gtexBrain, file = 'data/Reference/GTEx/GTEx_Brain_clinical.RDS')
saveRDS(gtexData, file = 'data/Reference/GTEx/GTEx_Brain_TPM.RDS')
