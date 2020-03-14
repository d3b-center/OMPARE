#####################
# Cleanup data for P5
#####################

# this code is to be run as is 
clinData.full <- clinData %>%
  filter(experimental_strategy == "RNA-Seq") %>%
  dplyr::select(Kids_First_Biospecimen_ID, sample_id, reported_gender, race, age_at_diagnosis_days, ethnicity, pathology_diagnosis, integrated_diagnosis, short_histology, broad_histology, primary_site)

clinData <- clinData %>%
  filter(experimental_strategy == "RNA-Seq") %>%
  dplyr::select(Kids_First_Biospecimen_ID, sample_id, pathology_diagnosis, integrated_diagnosis, short_histology, broad_histology, primary_site)
patVector <- c(sampleInfo$subjectID, sampleInfo$subjectID, rep(sampleInfo$tumorType, 4), sampleInfo$tumorLocation)
clinData <- rbind(clinData, patVector)
rownames(clinData) <- clinData$Kids_First_Biospecimen_ID

# Combine PBTA and Patient
combGenes <- intersect(rownames(res), rownames(expData))
res <- cbind(res[combGenes,], expData[combGenes,sampleInfo$subjectID])
colnames(res)[ncol(res)] <- sampleInfo$subjectID
resAll <- res

# Now remove genes that have less 20 FPKM
maxVals <- apply(res, FUN = max, MARGIN = 1)
res <- res[maxVals>50,]

# Order samples for expression and clinical file
common.smps <- intersect(colnames(res), rownames(clinData))
res <- res[,common.smps]
clinData <- clinData[common.smps,]

###########################
# Get Annotation data ready 
# Constrain columns
##########################

# Get most variable genes
# Top 10000
myCV <- function(x) { sd(x)/mean(x)}
myCVs <- apply(res, FUN=myCV, MARGIN=1)
resTmp <- res
resTmp["CV"] <- myCVs
resTmp <- resTmp[order(-resTmp[,"CV"]),]
resTmp <- resTmp[1:10000,]
resTmp <- resTmp[-ncol(resTmp)] # Remove cv

diseasetypes <- c("High-grade glioma", sampleInfo$tumorType)
clinDataHGG <- clinData[grepl(paste(diseasetypes, collapse = "|"), clinData$integrated_diagnosis),]
resTmpHGG <- resTmp[,rownames(clinDataHGG)]
allCor <- cor(resTmpHGG[ncol(resTmpHGG)], resTmpHGG)
allCor <- data.frame(t(allCor), check.names = F)
allCor[,"samps"] <- rownames(allCor)
allCor <- allCor[!grepl(sampleInfo$subjectID, rownames(allCor)),]
allCor <- allCor[intersect(rownames(allCor), survData$Kids_First_Biospecimen_ID),]
allCor <- allCor[order(-allCor[,1]),]


