# Script to create TCGA TPM matrix of unique gene symbols and colnames
library(RDiseaseXpress)
library(tidyverse)

# clinical file (n = 155)
tcga <- getSamples(myStudy = "TCGA")
tcga <- tcga[grep('Primary', tcga$definition),]
tcga.clinData <- tcga %>%
  filter(disease %in% "GBM") %>%
  mutate(short_histology = disease,
         broad_histology = disease_name,
         library_name = "polyA") %>%
  dplyr::select(sample_barcode, sample_id, gender, age_at_diagnosis_in_days, ethnicity, short_histology, broad_histology, overall_survival_time_in_days, vital_status, study_id, library_name) %>%
  mutate(vital_status = ifelse(vital_status == "dead", 1, 0)) %>%
  arrange(sample_id)

# add survival from pedcbio
tcga.gbm.pedcbio <-  read.delim('~/Projects/toil-rnaseq-20k/data/clinical/gbm_tcga_pub2013_clinical_data.tsv')
tcga.gbm.pedcbio <- tcga.gbm.pedcbio[,c('Patient.ID', 'Overall.Survival..Months.')]
tcga.gbm.clinData <- tcga.clinData
tcga.gbm.clinData$patient_barcode <- gsub('-[0-9]{2}A.*','',tcga.gbm.clinData$sample_barcode)
tcga.gbm.clinData <- merge(tcga.gbm.clinData, tcga.gbm.pedcbio, by.x = 'patient_barcode', by.y = 'Patient.ID', all.x = TRUE)
tcga.gbm.clinData$overall_survival_time_in_days <- round(tcga.gbm.clinData$Overall.Survival..Months.*30.417)
tcga.gbm.clinData$patient_barcode <- NULL
tcga.gbm.clinData$Overall.Survival..Months. <- NULL
tcga.gbm.clinData$overall_survival_time_in_days[is.na(tcga.gbm.clinData$overall_survival_time_in_days)] <- 'unavailable'
tcga.gbm.clinData <- tcga.gbm.clinData[order(tcga.gbm.clinData$sample_id),]

# TPM matrix
tcga.tpm <-  readRDS('data/Reference/TCGA/TCGA_matrix_TPM.RDS')
tcga.gbm.tpm <- tcga.tpm[,colnames(tcga.tpm) %in% tcga.gbm.clinData$sample_id]
tcga.gbm.tpm <- tcga.gbm.tpm[,order(colnames(tcga.gbm.tpm))]
tcga.gbm.clinData$sample_barcode <- gsub('R-.*','',tcga.gbm.clinData$sample_barcode)
if(identical(colnames(tcga.gbm.tpm), tcga.gbm.clinData$sample_id)){
  print("Proceed")
  colnames(tcga.gbm.tpm) <- tcga.gbm.clinData$sample_barcode
}

# save
system('mkdir -p data/Reference/TCGA/')
saveRDS(tcga.gbm.tpm, file = 'data/Reference/TCGA/TCGA_GBM_matrix_TPM.RDS')
saveRDS(tcga.gbm.clinData, file = "data/Reference/TCGA/TCGA_GBM_clinData.RDS")

