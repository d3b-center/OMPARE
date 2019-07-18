#####################
# Cleanup data for P5
#####################

# this code is to be run as is (will fix later)
clinData <- merge(clinData, mapping, by.x="Sample.ID", by.y="V1")
clinDataOrig <- clinData
clinData <- clinData[!grepl("CL", clinData[,"CBTTC_PAIRED_IDS"]),]
clinData <- clinData[,c("V2", "Cancer.Type", "Cancer.Type.Detailed", "TUMOR_TISSUE_SITE")]
clinData <- unique(clinData)
clinData <- clinData[!duplicated(clinData[,"V2"]),]
rownames(clinData)<- clinData[,"V2"]
patVector <- c("PatSample", "High-grade glioma/astrocytoma (WHO grade III/IV)", "High-grade glioma/astrocytoma (WHO grade III/IV)", "Unknown")
clinData <- rbind(clinData, patVector)
rownames(clinData)[nrow(clinData)] <- "PatSample"


# Filter data and get it by gene, remove all genes with "-", ".", "_"
rownames(res) <- res[,1]
combGenes <- intersect(rownames(res), rownames(expData))
res <- cbind(res[combGenes,], expData[combGenes,"FPKM"])
colnames(res)[ncol(res)] <- "PatSample"
res[,2] <- as.character(res[,2])
res <- res[!grepl("-", res[,2]),]
res <- res[!grepl("\\.", res[,2]),]
res <- res[!grepl("_", res[,2]),]
resAll <- res

# Now remove genes that have less 20 FPKM
maxVals <- apply(res[3:ncol(res)], FUN=max, MARGIN=1)
res <- res[maxVals>50,]

# Now take gene with max value
res[,"max"] <- apply(res[3:ncol(res)], FUN=max, MARGIN=1)
res <- res[order(-res[,"max"]),]
res <- res[!duplicated(res[,2]),]
rownames(res) <- res[,2]
res <- res[-1:-2]
res <- res[-ncol(res)] #Remove max
res <- res[,rownames(clinData)]

###########################
# Get Annotation data ready 
# Constrain columns
##########################

# Get most variable genes
myCV <- function(x) { sd(x)/mean(x)}
myCVs <- apply(res, FUN=myCV, MARGIN=1)
resTmp <- res
resTmp["CV"] <- myCVs
resTmp <- resTmp[order(-resTmp[,"CV"]),]
resTmp <- resTmp[1:10000,]
resTmp <- resTmp[-ncol(resTmp)] # Remove cv

clinDataHGG <- clinData[grepl("High-grade glioma", clinData[,2]),]
resTmpHGG <- resTmp[,rownames(clinDataHGG)]
allCor <- cor(resTmpHGG[ncol(resTmpHGG)], resTmpHGG)
allCor <- data.frame(t(allCor))
allCor[,"samps"] <- rownames(allCor)
allCor <- allCor[!grepl("Pat", rownames(allCor)),]
allCor <- allCor[intersect(rownames(allCor), survData[,"samps"]),]
allCor <- allCor[order(-allCor[,1]),]


