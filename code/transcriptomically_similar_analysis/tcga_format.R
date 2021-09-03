# script to format tcga + pnoc008 data

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
tcga_dir <- file.path(data_dir, "tcga")
module_dir <- file.path(root_dir, "code", "transcriptomically_similar_analysis")
output_dir <- file.path(patient_dir, "output", "transcriptomically_similar_analysis")
dir.create(output_dir, recursive = T, showWarnings = F)

# source functions
source(file.path(module_dir, "utils", "combine_and_batch_correct.R"))
source(file.path(module_dir, "utils", "get_data_for_immune_profile.R"))
source(file.path(module_dir, "utils", "get_most_variable_for_umap.R"))
source(file.path(module_dir, "utils", "extract_umap_nearest_neighbor_table.R"))
source(file.path(module_dir, "utils", "get_umap_output.R"))

# reference
# all PNOC008 patients (TPM matrix + clinical)
pnoc008_tpm <- readRDS(file.path(data_dir, 'pnoc008', 'pnoc008_tpm_matrix.rds'))
pnoc008_clinical <- readRDS(file.path(data_dir, 'pnoc008', 'pnoc008_clinical.rds'))

# TCGA GBM specific data (TPM)
tcga_gbm_tpm <- readRDS(file.path(data_dir, 'tcga', 'tcga_gbm_tpm_matrix.rds'))
tcga_gbm_clinical <- readRDS(file.path(data_dir, 'tcga', 'tcga_gbm_clinical.rds'))

# combine and correct tcga gbm + pnoc008
fname <- file.path(tcga_dir, "tcga_gbm_pnoc008_corrected_matrix.rds")
tcga_gbm_clinical <- tcga_gbm_clinical %>%
  mutate(pathology_diagnosis = short_histology,
         integrated_diagnosis = broad_histology)
res <- combine_and_batch_correct(expr1 = tcga_gbm_tpm, 
                                 clinical1 = tcga_gbm_clinical, 
                                 expr2 = pnoc008_tpm,
                                 clinical2 = pnoc008_clinical)
saveRDS(res, file = fname)
tcga_gbm_pnoc008_expr_uncorrected <- res$expr_uc
tcga_gbm_pnoc008_expr_corrected <- res$expr
tcga_gbm_pnoc008_clinical <- res$clinical

# get data for immune profile
tcga_gbm_pnoc008_immune_profile <- get_data_for_immune_profile(expr_uncorrected = tcga_gbm_pnoc008_expr_uncorrected, 
                                                               expr1 = tcga_gbm_tpm, 
                                                               patient_of_interest = patient)

# get 10000 most variable genes for umap plotting & nearest neighbor analysis
tcga_gbm_pnoc008_most_var <- get_most_variable_for_umap(expr_corrected = tcga_gbm_pnoc008_expr_corrected)

# nearest neighbor info using umap correlation
fname <- file.path(output_dir, "tcga_pnoc008_umap_output.rds")
tcga_gbm_pnoc008_umap_output <- get_umap_output(expr_most_var = tcga_gbm_pnoc008_most_var)  
saveRDS(tcga_gbm_pnoc008_umap_output, file = fname)

# embeddings: required for umap clustering plot
tcga_gbm_pnoc008_embedding <- as.data.frame(tcga_gbm_pnoc008_umap_output$embedding)

# extract nearest neighbor info
tcga_gbm_pnoc008_nn_table <- extract_umap_nearest_neighbor_table(umap_out = tcga_gbm_pnoc008_umap_output, 
                                                                 expr_most_var = tcga_gbm_pnoc008_most_var, 
                                                                 patient_of_interest = patient)

# mutational_analysis
tcga_gbm_pnoc008_nn_tpm <- tcga_gbm_pnoc008_expr_corrected[,colnames(tcga_gbm_pnoc008_expr_corrected) %in% tcga_gbm_pnoc008_nn_table$nearest_neighbor]

# required for pathway_analysis, kaplan meier, transcriptomically_similar analyses
fname <-  file.path(output_dir, "tcga_gbm_pnoc008_nn_table.rds")
tcga_gbm_pnoc008_nn_table <- tcga_gbm_pnoc008_nn_table[grep(patient, tcga_gbm_pnoc008_nn_table$nearest_neighbor, invert = TRUE),]
saveRDS(tcga_gbm_pnoc008_nn_table, file = fname)

