#####################
# Format TCGA data
#####################

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))
tcga_dir <- file.path(ref_dir, 'tcga')

# source functions
source(file.path(patient_level_analyses_utils, 'quiet.R'))
source(file.path(patient_level_analyses_utils, 'batch_correct.R'))

# TCGA specific
tcga_gbm_clinical <- tcga_gbm_clinical %>%
  dplyr::select(-c(overall_survival_time_in_days, vital_status)) %>%
  as.data.frame()
pnoc008_clinical_sub <- pnoc008_clinical[,c(rep('subjectID', 2), 'sex', 'age_diagnosis_days', 'ethnicity', rep('tumorType',2), 'study_id', 'library_name')] 
colnames(pnoc008_clinical_sub) <- colnames(tcga_gbm_clinical)
tcga_gbm_clinical <- rbind(tcga_gbm_clinical, pnoc008_clinical_sub)
rownames(tcga_gbm_clinical) <- tcga_gbm_clinical$sample_barcode

# Combine tcga.gbm and PNOC Patients
combGenes <- intersect(rownames(tcga_gbm_tpm), rownames(pnoc008_tpm))
tcga_gbm_tpm <- cbind(tcga_gbm_tpm[combGenes,], pnoc008_tpm[combGenes,])
tcga_gbm_tpm <- tcga_gbm_tpm[,rownames(tcga_gbm_clinical)]

# Correct for batch effect: study_id + library_name
tcga_gbm_clinical$batch <- paste0(tcga_gbm_clinical$study_id,'_', tcga_gbm_clinical$library_name)
fname <- file.path(tcga_dir, 'tcga_gbm_pnoc008_corrected_matrix.rds')
if(snv_pattern != "lancet" & file.exists(fname)){
  tcga_gbm_tpm <- readRDS(fname)
} else {
  tcga_gbm_tpm <- quiet(batch.correct(mat = tcga_gbm_tpm, clin = tcga_gbm_clinical))
  saveRDS(tcga_gbm_tpm, file = fname)
}

# keep full matrix for ImmuneProfile.R (only TCGA + PNOC patient of interest)
tcga_gbm_tpm_full <- tcga_gbm_tpm 
smps <- grep('TCGA', colnames(tcga_gbm_tpm_full), value = T)
smps <- c(smps, sampleInfo$subjectID)
tcga_gbm_tpm_all <- tcga_gbm_tpm_full[,colnames(tcga_gbm_tpm_full) %in% smps]

# Now remove genes that have max value < 20 TPM
maxVals <- apply(tcga_gbm_tpm, FUN = max, MARGIN = 1)
tcga_gbm_tpm <- tcga_gbm_tpm[maxVals>20,]

# Order samples for expression and clinical file
common.smps <- intersect(colnames(tcga_gbm_tpm), rownames(tcga_gbm_clinical))
tcga_gbm_tpm <- tcga_gbm_tpm[,common.smps]
tcga_gbm_clinical <- tcga_gbm_clinical[common.smps,]

# for dimensionality reduction visualization (dim_reduction_plot.R)
# Get top 1000 most variable genes
myCV <- function(x) { sd(x)/mean(x)}
myCVs <- apply(tcga_gbm_tpm, FUN=myCV, MARGIN=1)
tcga_gbm_tpm_tsne <- as.data.frame(tcga_gbm_tpm)
tcga_gbm_tpm_tsne$CV <- myCVs
tcga_gbm_tpm_tsne <- tcga_gbm_tpm_tsne[order(tcga_gbm_tpm_tsne$CV, decreasing = TRUE),]
if(nrow(tcga_gbm_tpm_tsne) >= 1000){
  tcga_gbm_tpm_tsne <- tcga_gbm_tpm_tsne[1:1000,]
}
tcga_gbm_tpm_tsne$CV <- NULL # Remove cv

# clustering using umap correlation
tcga_umap_output <- file.path(topDir, 'output', 'tcga_pnoc008_umap_output.rds')
if(file.exists(tcga_umap_output)){
  tcga_umap <- readRDS(file = tcga_umap_output)
} else {
  set.seed(100)
  tcga_umap <- uwot::umap(X = t(log2(tcga_gbm_tpm_tsne+1)), n_neighbors = 21, n_components = 2, metric = "correlation", ret_nn = TRUE, n_sgd_threads = 123L)
  
  # add colnames/rownames to embeddings
  colnames(tcga_umap$embedding) <- c("UMAP1", "UMAP2")
  rownames(tcga_umap$embedding) <- colnames(tcga_gbm_tpm_tsne)
  
  # add rownames to nearest neighbor
  rownames(tcga_umap$nn$correlation$idx) <- colnames(tcga_gbm_tpm_tsne)
  rownames(tcga_umap$nn$correlation$dist) <- colnames(tcga_gbm_tpm_tsne)
  
  # save output  
  saveRDS(tcga_umap, file = tcga_umap_output)
}

# embeddings: required for umap clustering plot
tcga_gbm_embedding <- as.data.frame(tcga_umap$embedding)

# extract nearest neighbor info
corr <- as.data.frame(tcga_umap$nn$correlation$idx) # nn
dist <- as.data.frame(tcga_umap$nn$correlation$dist) # distances
corr <- t(apply(corr, MARGIN = 1, FUN = function(x) colnames(tcga_gbm_tpm_tsne)[x]))
tcga_nn_table <- data.frame(nearest_neighbor = as.character(corr[grep(sampleInfo$subjectID, rownames(corr)),]), 
                            distance = as.numeric(dist[grep(sampleInfo$subjectID, rownames(dist)),]))
tcga_nn_table$distance <- round(tcga_nn_table$distance, digits = 3)

# required for pathway_analysis, kaplan meier, transcriptomically_similar analyses
tcga_gbm_allcor <- tcga_nn_table[grep(sampleInfo$subjectID, tcga_nn_table$nearest_neighbor, invert = TRUE),]

# immune_profile, ssgsea, mutational_analysis
# tcga_gbm_topcor <- tcga_gbm_tpm_full[,colnames(tcga_gbm_tpm_full) %in% tcga_nn_table$nearest_neighbor]
