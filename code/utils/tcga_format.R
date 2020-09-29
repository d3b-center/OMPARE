#####################
# Format TCGA data
#####################

# TCGA specific
tcga.gbm.clinData <- tcga.gbm.clinData %>%
  dplyr::select(-c(overall_survival_time_in_days, vital_status)) %>%
  as.data.frame()
pat.clinData <- pnoc008.clinData[,c(rep('subjectID', 2), 'sex', 'age_diagnosis_days', 'ethnicity', rep('tumorType',2), 'study_id', 'library_name')] 
colnames(pat.clinData) <- colnames(tcga.gbm.clinData)
tcga.gbm.clinData <- rbind(tcga.gbm.clinData, pat.clinData)
rownames(tcga.gbm.clinData) <- tcga.gbm.clinData$sample_barcode

# Combine tcga.gbm and PNOC Patients
combGenes <- intersect(rownames(tcga.gbm.mat), rownames(pnoc008.data))
tcga.gbm.mat <- cbind(tcga.gbm.mat[combGenes,], pnoc008.data[combGenes,])
tcga.gbm.mat <- tcga.gbm.mat[,rownames(tcga.gbm.clinData)]

# Correct for batch effect: study_id + library_name
tcga.gbm.clinData$batch <- paste0(tcga.gbm.clinData$study_id,'_', tcga.gbm.clinData$library_name)
if(snv_pattern != "lancet" & file.exists('data/Reference/TCGA/tcga_gbm_pnoc008_corrected_matrix.rds')){
  tcga.gbm.mat <- readRDS('data/Reference/TCGA/tcga_gbm_pnoc008_corrected_matrix.rds')
} else {
  tcga.gbm.mat <- quiet(batch.correct(mat = tcga.gbm.mat, clin = tcga.gbm.clinData))
  saveRDS(tcga.gbm.mat, file = 'data/Reference/TCGA/tcga_gbm_pnoc008_corrected_matrix.rds')
}

# keep full matrix for ImmuneProfile.R (only TCGA + PNOC patient of interest)
tcga.gbm.mat.full <- tcga.gbm.mat 
smps <- grep('TCGA_', colnames(tcga.gbm.mat.full), value = T)
smps <- c(smps, sampleInfo$subjectID)
tcga.gbm.mat.all <- tcga.gbm.mat.full[,colnames(tcga.gbm.mat.full) %in% smps]

# Now remove genes that have max value < 20 TPM
maxVals <- apply(tcga.gbm.mat, FUN = max, MARGIN = 1)
tcga.gbm.mat <- tcga.gbm.mat[maxVals>20,]

# Order samples for expression and clinical file
common.smps <- intersect(colnames(tcga.gbm.mat), rownames(tcga.gbm.clinData))
tcga.gbm.mat <- tcga.gbm.mat[,common.smps]
tcga.gbm.clinData <- tcga.gbm.clinData[common.smps,]

###########################
# Get Annotation data ready 
# Constrain columns
##########################

# for getTSNEPlot.R
# Get top 10000 most variable genes
myCV <- function(x) { sd(x)/mean(x)}
myCVs <- apply(tcga.gbm.mat, FUN=myCV, MARGIN=1)
tcga.gbm.mat.tsne <- as.data.frame(tcga.gbm.mat)
tcga.gbm.mat.tsne$CV <- myCVs
tcga.gbm.mat.tsne <- tcga.gbm.mat.tsne[order(tcga.gbm.mat.tsne$CV, decreasing = TRUE),]
if(nrow(tcga.gbm.mat.tsne) >= 10000){
  tcga.gbm.mat.tsne <- tcga.gbm.mat.tsne[1:10000,]
}
tcga.gbm.mat.tsne$CV <- NULL # Remove cv

# for clustering
# use UMAP correlation
tcga.umap.output <- file.path(topDir, 'Summary/tcga_pnoc008_umap_output.rds')
if(file.exists(tcga.umap.output)){
  tcga.umap <- readRDS(file = tcga.umap.output)
} else {
  set.seed(100)
  tcga.umap <- uwot::umap(X = t(log2(tcga.gbm.mat.tsne+1)), n_neighbors = 21, n_components = 2, metric = "correlation", ret_nn = TRUE, n_sgd_threads = 123L)
  saveRDS(tcga.umap, file = tcga.umap.output)
}
tcga.gbm.embedding <- as.data.frame(tcga.umap$embedding)
colnames(tcga.gbm.embedding) <- c("UMAP1", "UMAP2")

# for getKMPlot.R and getSimilarPatients.R
# extract nearest neighbor info
corr <- as.data.frame(tcga.umap$nn$correlation$idx) # nn
dist <- as.data.frame(tcga.umap$nn$correlation$dist) # distances
corr <- t(apply(corr, MARGIN = 1, FUN = function(x) colnames(tcga.gbm.mat.tsne)[x]))
rownames(corr) <- colnames(tcga.gbm.mat.tsne)
rownames(dist) <- colnames(tcga.gbm.mat.tsne)
tcga_nn_table <- data.frame(nearest_neighbor = as.character(corr[grep(sampleInfo$subjectID, rownames(corr)),]), 
                       distance = as.numeric(dist[grep(sampleInfo$subjectID, rownames(dist)),]))
tcga_nn_table$distance <- round(tcga_nn_table$distance, digits = 3)
tcga.gbm.allCor <- tcga_nn_table[grep(sampleInfo$subjectID, tcga_nn_table$nearest_neighbor, invert = TRUE),]

# Immune profile, ssGSEA, recurrent alterations (keep POI)
tcga.gbm.topCor <- tcga.gbm.mat.full[,colnames(tcga.gbm.mat.full) %in% tcga_nn_table$nearest_neighbor]
