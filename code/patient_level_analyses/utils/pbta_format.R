#####################
# Format PBTA data
#####################

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))
pbta_dir <- file.path(ref_dir, 'pbta')

# source functions
source(file.path(patient_level_analyses_utils, 'quiet.R'))
source(file.path(patient_level_analyses_utils, 'batch_correct.R'))

# PBTA clinical
pbta_clinical <- pbta_clinical %>%
  filter(experimental_strategy == "RNA-Seq") %>%
  mutate(sample_barcode = Kids_First_Biospecimen_ID,
         study_id = "PBTA",
         library_name = RNA_library) %>%
  dplyr::select(sample_barcode, sample_id, reported_gender, age_at_diagnosis_days, ethnicity, pathology_diagnosis, integrated_diagnosis, short_histology, broad_histology, primary_site, study_id, library_name)

# PNOC008 clinical 
# merge with PBTA clinical
pnoc008_clinical_sub <- pnoc008_clinical[,c(rep('subjectID', 2), 'sex', 'age_diagnosis_days', 'ethnicity', rep('tumorType', 4), 'tumorLocation', 'study_id', 'library_name')] 
colnames(pnoc008_clinical_sub) <- colnames(pbta_clinical)
pbta_clinical <- rbind(pbta_clinical, pnoc008_clinical_sub)
rownames(pbta_clinical) <- pbta_clinical$sample_barcode

# Combine PBTA and PNOC Patients expression matrix
combGenes <- intersect(rownames(pbta_tpm), rownames(pnoc008_tpm))
pbta_tpm <- cbind(pbta_tpm[combGenes,], pnoc008_tpm[combGenes,])
pbta_tpm <- pbta_tpm[,rownames(pbta_clinical)]

# Correct for batch effect: study_id + library_name
pbta_clinical$batch <- paste0(pbta_clinical$study_id,'_', pbta_clinical$library_name)
fname <- file.path(pbta_dir, 'pbta_pnoc008_corrected_matrix.rds')
if(snv_pattern != "lancet" & file.exists(fname)){
  pbta_tpm <- readRDS(fname)
} else {
  pbta_tpm <- quiet(batch.correct(mat = pbta_tpm, clin = pbta_clinical))
  saveRDS(pbta_tpm, file = fname)
}

# keep full matrix for ImmuneProfile.R (only PBTA + PNOC patient of interest)
pbta_tpm_full <- pbta_tpm 
smps <- grep('BS_', colnames(pbta_tpm_full), value = T)
smps <- c(smps, sampleInfo$subjectID)
pbta_tpm_all <- pbta_tpm_full[,colnames(pbta_tpm_full) %in% smps]

# Now remove genes that have max value < 20 TPM
maxVals <- apply(pbta_tpm, FUN = max, MARGIN = 1)
pbta_tpm <- pbta_tpm[maxVals > 20,]

# Order samples for expression and clinical file
common.smps <- intersect(colnames(pbta_tpm), rownames(pbta_clinical))
pbta_tpm <- pbta_tpm[,common.smps]
pbta_clinical <- pbta_clinical[common.smps,]

# for dimensionality reduction visualization (dim_reduction_plot.R)
# Get top 1000 most variable genes
myCV <- function(x) { sd(x)/mean(x)}
myCVs <- apply(pbta_tpm, FUN = myCV, MARGIN=1)
pbta_tpm_tsne <- as.data.frame(pbta_tpm)
pbta_tpm_tsne$CV <- myCVs
pbta_tpm_tsne <- pbta_tpm_tsne[order(pbta_tpm_tsne$CV, decreasing = TRUE),]
if(nrow(pbta_tpm_tsne) >= 1000){
  pbta_tpm_tsne <- pbta_tpm_tsne[1:1000,]
}
pbta_tpm_tsne$CV <- NULL # Remove cv

# clustering using umap correlation
pbta_umap_output <- file.path(topDir, 'output', 'pbta_pnoc008_umap_output.rds')
if(file.exists(pbta_umap_output)){
  pbta_umap <- readRDS(file = pbta_umap_output)
} else {
  set.seed(100)
  pbta_umap <- uwot::umap(X = t(log2(pbta_tpm_tsne+1)), n_neighbors = 21, n_components = 2, metric = "correlation", ret_nn = TRUE, n_sgd_threads = 123L)
  
  # add colnames/rownames to embeddings
  colnames(pbta_umap$embedding) <- c("UMAP1", "UMAP2")
  rownames(pbta_umap$embedding) <- colnames(pbta_tpm_tsne)
  
  # add rownames to nearest neighbor
  rownames(pbta_umap$nn$correlation$idx) <- colnames(pbta_tpm_tsne)
  rownames(pbta_umap$nn$correlation$dist) <- colnames(pbta_tpm_tsne)
  
  # save output
  saveRDS(pbta_umap, file = pbta_umap_output)
}

# embeddings: required for umap clustering plot
pbta_embedding <- as.data.frame(pbta_umap$embedding)

# extract nearest neighbor info
corr <- as.data.frame(pbta_umap$nn$correlation$idx) # nn
dist <- as.data.frame(pbta_umap$nn$correlation$dist) # distances
corr <- t(apply(corr, MARGIN = 1, FUN = function(x) colnames(pbta_tpm_tsne)[x]))
pbta_nn_table <- data.frame(nearest_neighbor = as.character(corr[grep(sampleInfo$subjectID, rownames(corr)),]), 
                            distance = as.numeric(dist[grep(sampleInfo$subjectID, rownames(dist)),]))
pbta_nn_table$distance <- round(pbta_nn_table$distance, digits = 3)

# required for pathway_analysis, kaplan meier, transcriptomically_similar 
pbta_allcor <- pbta_nn_table[grep(sampleInfo$subjectID, pbta_nn_table$nearest_neighbor, invert = TRUE),]

# HGAT samples + POI (ssGSEA)
# pbta_hgg_clinical <- pbta_clinical %>% filter(short_histology == "HGAT")
# pbta_hgg_clinical <- c(sampleInfo$subjectID, pbta_hgg_clinical$sample_barcode)
# pbta_hgg_tpm <- pbta_tpm_full[,colnames(pbta_tpm_full) %in% pbta_hgg_clinical]

# immune_profile_topcor, ssgsea, mutational_analysis
pbta_topcor <- pbta_tpm_full[,colnames(pbta_tpm_full) %in% pbta_nn_table$nearest_neighbor]
