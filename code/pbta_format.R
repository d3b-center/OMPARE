#####################
# Format PBTA data
#####################

# PBTA specific
# pbta.clinData.full <- pbta.clinData %>%
pbta.clinData <- pbta.clinData %>%
  filter(experimental_strategy == "RNA-Seq") %>%
  mutate(sample_barcode = Kids_First_Biospecimen_ID,
         study_id = "PBTA") %>%
  dplyr::select(sample_barcode, sample_id, reported_gender, age_at_diagnosis_days, ethnicity, pathology_diagnosis, integrated_diagnosis, short_histology, broad_histology, primary_site, study_id)

pat.clinData <- pnoc008.clinData[,c(rep('subjectID', 2), 'sex', 'AgeAtCollection', 'ethnicity', rep('tumorType', 4), 'tumorLocation', 'study_id')] 
colnames(pat.clinData) <- colnames(pbta.clinData)
pbta.clinData <- rbind(pbta.clinData, pat.clinData)
rownames(pbta.clinData) <- pbta.clinData$sample_barcode

# Combine PBTA and PNOC Patients
combGenes <- intersect(rownames(pbta.mat), rownames(pnoc008.data))
pbta.mat <- cbind(pbta.mat[combGenes,], pnoc008.data[combGenes,])

# keep full matrix for ImmuneProfile.R (only PBTA + PNOC patient of interest)
pbta.mat.full <- pbta.mat 
smps <- grep('BS_', colnames(pbta.mat.full), value = T)
smps <- c(smps, sampleInfo$subjectID)
pbta.mat.all <- pbta.mat.full[,colnames(pbta.mat.full) %in% smps]

# Now remove genes that have max value < 50 TPM
maxVals <- apply(pbta.mat, FUN = max, MARGIN = 1)
pbta.mat <- pbta.mat[maxVals>50,]

# Order samples for expression and clinical file
common.smps <- intersect(colnames(pbta.mat), rownames(pbta.clinData))
pbta.mat <- pbta.mat[,common.smps]
pbta.clinData <- pbta.clinData[common.smps,]

###########################
# Get Annotation data ready 
# Constrain columns
##########################

# for getTSNEPlot.R
# Get top 10000 most variable genes
myCV <- function(x) { sd(x)/mean(x)}
myCVs <- apply(pbta.mat, FUN=myCV, MARGIN=1)
pbta.mat.tsne <- pbta.mat
pbta.mat.tsne["CV"] <- myCVs
pbta.mat.tsne <- pbta.mat.tsne[order(-pbta.mat.tsne[,"CV"]),]
pbta.mat.tsne <- pbta.mat.tsne[1:10000,]
pbta.mat.tsne <- pbta.mat.tsne[-ncol(pbta.mat.tsne)] # Remove cv

# for getKMPlot.R and getSimilarPatients.R
pbta.allCor <- cor(pbta.mat.tsne[sampleInfo$subjectID], pbta.mat.tsne)
pbta.allCor <- data.frame(t(pbta.allCor), check.names = F)
pbta.allCor[,"sample_barcode"] <- rownames(pbta.allCor)
pbta.allCor <- pbta.allCor[!grepl(sampleInfo$subjectID, rownames(pbta.allCor)),]
pbta.allCor <- pbta.allCor[order(-pbta.allCor[,1]),]
pbta.allCor[,1] <- round(pbta.allCor[,1], 3)

# get matrix of top 20 correlated samples (for Immune profile of genomically similar patients)
pbta.topCor <- pbta.allCor[1:20,'sample_barcode']
pbta.topCor <- c(pbta.topCor, sampleInfo$subjectID) # add patient of interest
pbta.topCor <- pbta.mat.full[,colnames(pbta.mat.full) %in% pbta.topCor]
