# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))
source(file.path(patient_level_analyses_utils, "data_formatting_functions.R"))
pbta_dir <- file.path(ref_dir, 'pbta')

# combine and correct pbta + pnoc008
fname <- file.path(pbta_dir, "pbta_pnoc008_corrected_matrix.rds")
if(file.exists(fname) & snv_pattern != "lancet"){
  res <- readRDS(fname)
} else {
  pbta_clinical <- pbta_clinical %>%
    filter(experimental_strategy == "RNA-Seq") %>%
    mutate(sample_barcode = Kids_First_Biospecimen_ID,
           study_id = "PBTA",
           library_name = RNA_library,
           gender = reported_gender, 
           age_at_diagnosis_in_days = age_at_diagnosis_days)
  res <- combine_and_batch_correct(expr1 = pbta_tpm, expr2 = pnoc008_tpm,  clinical1 = pbta_clinical, clinical2 = pnoc008_clinical)
  saveRDS(res, file = fname)
}
pbta_pnoc008_expr_uncorrected <- res$expr_uc
pbta_pnoc008_expr_corrected <- res$expr
pbta_pnoc008_clinical <- res$clinical

# get data for immune profile
pbta_pnoc008_immune_profile <- get_data_for_immune_profile(expr_uncorrected = pbta_pnoc008_expr_uncorrected, expr1 = pbta_tpm, sampleInfo = sampleInfo)

# get 10000 most variable genes for umap plotting & nearest neighbor analysis
pbta_pnoc008_most_var <- get_most_variable_for_umap(expr_corrected = pbta_pnoc008_expr_corrected)

# nearest neighbor info using umap correlation
fname <- file.path(topDir, "output", "pbta_pnoc008_umap_output.rds")
if(file.exists(fname)){
  pbta_pnoc008_umap_output <- readRDS(fname)
} else {
  pbta_pnoc008_umap_output <- get_umap_output(expr_most_var = pbta_pnoc008_most_var)  
  saveRDS(pbta_pnoc008_umap_output, file = fname)
}

# embeddings: required for umap clustering plot
pbta_pnoc008_embedding <- as.data.frame(pbta_pnoc008_umap_output$embedding)

# extract nearest neighbor info
pbta_pnoc008_nn_table <- extract_umap_nearest_neighbor_table(umap_out = pbta_pnoc008_umap_output, expr_most_var = pbta_pnoc008_most_var, sampleInfo = sampleInfo)

# immune_profile_topcor, ssgsea, mutational_analysis
pbta_pnoc008_nn_tpm <- pbta_pnoc008_expr_corrected[,colnames(pbta_pnoc008_expr_corrected) %in% pbta_pnoc008_nn_table$nearest_neighbor]

# required for pathway_analysis, kaplan meier, transcriptomically_similar analyses
pbta_pnoc008_nn_table <- pbta_pnoc008_nn_table[grep(sampleInfo$subjectID, pbta_pnoc008_nn_table$nearest_neighbor, invert = TRUE),]
