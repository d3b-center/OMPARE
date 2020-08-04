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
tcga.gbm.mat <- quiet(batch.correct(mat = tcga.gbm.mat, clin = tcga.gbm.clinData))

# keep full matrix for ImmuneProfile.R (not required for TCGA for now)
# tcga.gbm.mat.all <- tcga.gbm.mat 

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
set.seed(100)
ump <- uwot::umap(X = t(log2(tcga.gbm.mat.tsne+1)), n_neighbors = 21, n_components = 2, metric = "correlation", ret_nn = TRUE, n_sgd_threads = 123L)
embedding <- as.data.frame(ump$embedding)
colnames(embedding) <- c("UMAP1", "UMAP2")
tcga.gbm.embedding <- embedding

# for getKMPlot.R and getSimilarPatients.R
# extract nearest neighbor info
corr <- as.data.frame(ump$nn$correlation$idx) # nn
dist <- as.data.frame(ump$nn$correlation$dist) # distances
corr <- t(apply(corr, MARGIN = 1, FUN = function(x) colnames(tcga.gbm.mat.tsne)[x]))
rownames(corr) <- colnames(tcga.gbm.mat.tsne)
rownames(dist) <- colnames(tcga.gbm.mat.tsne)
nn_table <- data.frame(nearest_neighbor = as.character(corr[grep(sampleInfo$subjectID, rownames(corr)),]), 
                       distance = as.numeric(dist[grep(sampleInfo$subjectID, rownames(dist)),]))
nn_table$distance <- round(nn_table$distance, digits = 3)
tcga.gbm.allCor <- nn_table[grep(sampleInfo$subjectID, nn_table$nearest_neighbor, invert = TRUE),]

# for getKMPlot.R and getSimilarPatients.R
# tcga.gbm.allCor <- cor(x = tcga.gbm.mat.tsne[sampleInfo$subjectID], y = tcga.gbm.mat.tsne)
# tcga.gbm.allCor <- data.frame(t(tcga.gbm.allCor), check.names = F)
# tcga.gbm.allCor[,"sample_barcode"] <- rownames(tcga.gbm.allCor)
# tcga.gbm.allCor <- tcga.gbm.allCor[!grepl(sampleInfo$subjectID, rownames(tcga.gbm.allCor)),]
# tcga.gbm.allCor <- tcga.gbm.allCor[order(tcga.gbm.allCor[,1], decreasing = TRUE),]
# tcga.gbm.allCor[,1] <- round(tcga.gbm.allCor[,1], 3)
