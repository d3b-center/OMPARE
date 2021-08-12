# script to format tcga data

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))
source(file.path(patient_level_analyses_utils, "data_formatting_functions.R"))
tcga_dir <- file.path(ref_dir, 'tcga')

# combine and correct tcga gbm + pnoc008
fname <- file.path(tcga_dir, "tcga_gbm_pnoc008_corrected_matrix.rds")
if(file.exists(fname) & snv_caller != "lancet"){
  res <- readRDS(fname)
} else {
  tcga_gbm_clinical <- tcga_gbm_clinical %>%
    mutate(pathology_diagnosis = short_histology,
           integrated_diagnosis = broad_histology)
  res <- combine_and_batch_correct(expr1 = tcga_gbm_tpm, expr2 = pnoc008_tpm, clinical1 = tcga_gbm_clinical, clinical2 = pnoc008_clinical)
  saveRDS(res, file = fname)
}
tcga_gbm_pnoc008_expr_uncorrected <- res$expr_uc
tcga_gbm_pnoc008_expr_corrected <- res$expr
tcga_gbm_pnoc008_clinical <- res$clinical

# get data for immune profile
tcga_gbm_pnoc008_immune_profile <- get_data_for_immune_profile(expr_uncorrected = tcga_gbm_pnoc008_expr_uncorrected, expr1 = tcga_gbm_tpm, sampleInfo = sampleInfo)

# get 10000 most variable genes for umap plotting & nearest neighbor analysis
tcga_gbm_pnoc008_most_var <- get_most_variable_for_umap(expr_corrected = tcga_gbm_pnoc008_expr_corrected)

# nearest neighbor info using umap correlation
fname <- file.path(topDir, "output", "tcga_pnoc008_umap_output.rds")
if(file.exists(fname)){
  tcga_gbm_pnoc008_umap_output <- readRDS(fname)
} else {
  tcga_gbm_pnoc008_umap_output <- get_umap_output(expr_most_var = tcga_gbm_pnoc008_most_var)  
  saveRDS(tcga_gbm_pnoc008_umap_output, file = fname)
}

# embeddings: required for umap clustering plot
tcga_gbm_pnoc008_embedding <- as.data.frame(tcga_gbm_pnoc008_umap_output$embedding)

# extract nearest neighbor info
tcga_gbm_pnoc008_nn_table <- extract_umap_nearest_neighbor_table(umap_out = tcga_gbm_pnoc008_umap_output, expr_most_var = tcga_gbm_pnoc008_most_var, sampleInfo = sampleInfo)

# mutational_analysis
tcga_gbm_pnoc008_nn_tpm <- tcga_gbm_pnoc008_expr_corrected[,colnames(tcga_gbm_pnoc008_expr_corrected) %in% tcga_gbm_pnoc008_nn_table$nearest_neighbor]

# required for pathway_analysis, kaplan meier, transcriptomically_similar analyses
tcga_gbm_pnoc008_nn_table <- tcga_gbm_pnoc008_nn_table[grep(sampleInfo$subjectID, tcga_gbm_pnoc008_nn_table$nearest_neighbor, invert = TRUE),]

