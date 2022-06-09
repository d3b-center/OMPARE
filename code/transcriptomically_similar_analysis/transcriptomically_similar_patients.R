# script to combine datasets and identify transcriptomically similar samples to patient of interest

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
module_dir <- file.path(root_dir, "code", "transcriptomically_similar_analysis")

# source functions
source(file.path(module_dir, "utils", "combine_and_batch_correct.R"))
source(file.path(module_dir, "utils", "get_data_for_immune_profile.R"))
source(file.path(module_dir, "utils", "get_most_variable_for_umap.R"))
source(file.path(module_dir, "utils", "extract_umap_nearest_neighbor_table.R"))
source(file.path(module_dir, "utils", "get_umap_output.R"))

# function
tns_similar_analysis <- function(ref_cancer_dir, patient_dir, sample_info, tpm_data, prefix){
  
  # output directory
  output_dir <- file.path(patient_dir, "output", "transcriptomically_similar_analysis")
  dir.create(output_dir, recursive = T, showWarnings = F)
  
  # patient_of_interest
  patient_of_interest <- sample_info %>%
    filter(experimental_strategy == "RNA-Seq") %>%
    .$Kids_First_Biospecimen_ID
  
  # reference histology
  ref_clinical <- list.files(path = ref_cancer_dir, pattern = "*.tsv", full.names = T)
  ref_clinical <- data.table::fread(ref_clinical)
  ref_clinical <- ref_clinical %>%
    filter(experimental_strategy == "RNA-Seq")
  
  # remove patient of interest
  if(length(intersect(ref_clinical$Kids_First_Biospecimen_ID, patient_of_interest)) != 0){
    ref_clinical <- ref_clinical %>%
      filter(!Kids_First_Biospecimen_ID %in% patient_of_interest)
  }
  
  # read tpm
  ref_tpm <- list.files(path = ref_cancer_dir, pattern = "*tpm.rds", full.names = T)
  ref_tpm <- readRDS(ref_tpm)
  
  # remove patient of interest
  if(length(intersect(colnames(ref_tpm), patient_of_interest)) != 0){
    ref_tpm <- ref_tpm %>%
      dplyr::select(-c(patient_of_interest))
  }
  
  # use common ids
  common_cols <- intersect(colnames(ref_tpm), ref_clinical$Kids_First_Biospecimen_ID)
  ref_clinical <- ref_clinical %>%
    filter(Kids_First_Biospecimen_ID %in% common_cols)
  ref_tpm <- ref_tpm[,ref_clinical$Kids_First_Biospecimen_ID]
  
  # combine and correct
  res <- combine_and_batch_correct(ref_expr = ref_tpm,  
                                   ref_clinical = ref_clinical, 
                                   subject_expr = tpm_data,
                                   subject_clinical = sample_info %>% filter(experimental_strategy == "RNA-Seq"))
  
  expr_uncorrected <- res$expr_uncorrected
  expr_corrected <- res$expr_corrected
  clinical <- res$clinical
  saveRDS(clinical, file = file.path(output_dir, paste0(prefix, '_patient_combined_clinical_input.rds')))
  
  # get data for immune profile
  immune_profile <- get_data_for_immune_profile(expr_uncorrected = expr_uncorrected, 
                                                ref_expr = ref_tpm, 
                                                patient_of_interest = patient_of_interest)
  saveRDS(immune_profile, file = file.path(output_dir, paste0(prefix, '_immune_profile_input.rds')))
  
  # get 10000 most variable genes for umap plotting & nearest neighbor analysis
  expr_most_var <- get_most_variable_for_umap(expr_corrected = expr_corrected)
  
  # nearest neighbor info using umap correlation
  fname <- file.path(output_dir, paste0(prefix, "_umap_output.rds"))
  umap_output <- get_umap_output(expr_most_var = expr_most_var)  
  saveRDS(umap_output, file = fname)
  
  # embeddings: required for umap clustering plot
  fname <- file.path(output_dir, paste0(prefix, "_umap_embedding.rds"))
  embedding <- as.data.frame(umap_output$embedding)
  saveRDS(embedding, file = fname)
  
  # extract nearest neighbor info
  nn_table <- extract_umap_nearest_neighbor_table(umap_out = umap_output, 
                                                  expr_most_var = expr_most_var, 
                                                  patient_of_interest = patient_of_interest)
  
  # immune_profile_topcor, ssgsea, mutational_analysis
  fname <-  file.path(output_dir, paste0(prefix, "_nn_tpm.rds"))
  nn_tpm <- expr_corrected[,colnames(expr_corrected) %in% nn_table$nearest_neighbor]
  saveRDS(nn_tpm, file = fname)
  
  # required for pathway_analysis, transcriptomically_similar analyses
  fname <-  file.path(output_dir, paste0(prefix, "_nn_table.rds"))
  nn_table <- nn_table[grep(patient_of_interest, nn_table$nearest_neighbor, invert = TRUE),]
  saveRDS(nn_table, file = fname)
}


