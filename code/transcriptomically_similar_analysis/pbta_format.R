# script to format pbta + pnoc008 data

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
pbta_dir <- file.path(data_dir, "pbta")
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

# PBTA specific mRNA data (TPM)
pbta_full_tpm <- readRDS(file.path(data_dir, 'pbta', 'pbta-gene-expression-rsem-tpm-collapsed.polya.stranded.rds'))
pbta_tpm <- pbta_full_tpm
pbta_clinical <- read.delim(file.path(data_dir, 'pbta', 'pbta-histologies.tsv'), stringsAsFactors = F)

# combine and correct pbta + pnoc008
fname <- file.path(pbta_dir, "pbta_pnoc008_corrected_matrix.rds")
pbta_clinical <- pbta_clinical %>%
  filter(experimental_strategy == "RNA-Seq",
         Kids_First_Biospecimen_ID %in% colnames(pbta_tpm)) %>% # will remove this after updating everything to v18
  mutate(study_id = "PBTA",
         library_name = RNA_library,
         gender = reported_gender, 
         age_at_diagnosis_in_days = age_at_diagnosis_days)
res <- combine_and_batch_correct(expr1 = pbta_tpm,  
                                 clinical1 = pbta_clinical, 
                                 expr2 = pnoc008_tpm,
                                 clinical2 = pnoc008_clinical)
saveRDS(res, file = fname)
pbta_pnoc008_expr_uncorrected <- res$expr_uc
pbta_pnoc008_expr_corrected <- res$expr
pbta_pnoc008_clinical <- res$clinical

# get data for immune profile
pbta_pnoc008_immune_profile <- get_data_for_immune_profile(expr_uncorrected = pbta_pnoc008_expr_uncorrected, 
                                                           expr1 = pbta_tpm, 
                                                           patient_of_interest = patient)

# get 10000 most variable genes for umap plotting & nearest neighbor analysis
pbta_pnoc008_most_var <- get_most_variable_for_umap(expr_corrected = pbta_pnoc008_expr_corrected)

# nearest neighbor info using umap correlation
fname <- file.path(output_dir, "pbta_pnoc008_umap_output.rds")
pbta_pnoc008_umap_output <- get_umap_output(expr_most_var = pbta_pnoc008_most_var)  
saveRDS(pbta_pnoc008_umap_output, file = fname)

# embeddings: required for umap clustering plot
pbta_pnoc008_embedding <- as.data.frame(pbta_pnoc008_umap_output$embedding)

# extract nearest neighbor info
pbta_pnoc008_nn_table <- extract_umap_nearest_neighbor_table(umap_out = pbta_pnoc008_umap_output, 
                                                             expr_most_var = pbta_pnoc008_most_var, 
                                                             patient_of_interest = patient)

# immune_profile_topcor, ssgsea, mutational_analysis
pbta_pnoc008_nn_tpm <- pbta_pnoc008_expr_corrected[,colnames(pbta_pnoc008_expr_corrected) %in% pbta_pnoc008_nn_table$nearest_neighbor]

# required for pathway_analysis, transcriptomically_similar analyses
fname <-  file.path(output_dir, "pbta_pnoc008_nn_table.rds")
pbta_pnoc008_nn_table <- pbta_pnoc008_nn_table[grep(patient, pbta_pnoc008_nn_table$nearest_neighbor, invert = TRUE),]
saveRDS(pbta_pnoc008_nn_table, file = fname)

# kaplan meier survival curves
# subset to top 20 HGAT only 
pbta_hgat_pnoc008_clinical <- pbta_pnoc008_clinical %>% 
  filter(short_histology == "HGAT",
         study_id == "PBTA")

# subset expression
pbta_hgat_pnoc008_expr_corrected <- pbta_pnoc008_expr_corrected %>%
  as.data.frame() %>%
  dplyr::select(patient, pbta_hgat_pnoc008_clinical$kids_first_biospecimen_id)

# get 10000 most variable genes for umap plotting & nearest neighbor analysis
pbta_hgat_pnoc008_most_var <- get_most_variable_for_umap(expr_corrected = pbta_hgat_pnoc008_expr_corrected)

# nearest neighbor info using umap correlation
fname <- file.path(output_dir, "pbta_hgat_pnoc008_umap_output.rds")
pbta_hgat_pnoc008_umap_output <- get_umap_output(expr_most_var = pbta_hgat_pnoc008_most_var)  
saveRDS(pbta_pnoc008_umap_output, file = fname)

# extract nearest neighbor info (top 20 nearest HGAT samples)
fname <- file.path(output_dir, "pbta_hgat_pnoc008_nn_table.rds")
pbta_hgat_pnoc008_nn_table <- extract_umap_nearest_neighbor_table(umap_out = pbta_hgat_pnoc008_umap_output, 
                                                                  expr_most_var = pbta_hgat_pnoc008_most_var, 
                                                                  patient_of_interest = patient)
pbta_hgat_pnoc008_nn_table <- pbta_hgat_pnoc008_nn_table[grep(patient, pbta_hgat_pnoc008_nn_table$nearest_neighbor, invert = TRUE),]
saveRDS(pbta_hgat_pnoc008_nn_table, file = fname)
