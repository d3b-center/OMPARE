##########################
# Functions required for :
# Key Clinical Findings P1
##########################


patientSampleInfo <- function() {
  df1 <- data.frame(c("Subject ID", "Sex", "Age", "Ethnicity"), c(getSubjectID(), getSex(), getAge(), getEthnicity()))
  df2 <- data.frame(c("Medical Facility", "Primary Physician", "Pathologist", "Lab Director"), c(getMedicalFacility(), getPrimPhysician(), getPathologist(), getLabDirector()))
  df3 <- data.frame(c("Collection Date", "Tumor Location", "Tumor Type", "P/R"), c(getCollectionDate(), getTumorLocation(), getTumorType(), getPrimRelapse()))
  return(cbind(df1,df2, df3))
}

keyClinicalFindingsTable <- function() {
  return(highConfidenceFindingsTable())
}

diseaseSpecificInformation <- function() {
  tmpGeneFindings <- allFindingsTable()
  if(nrow(tmpGeneFindings) > 0){
    tmpGeneFindings <- tmpGeneFindings[!grepl("Pathway", tmpGeneFindings[,"Type"]),]
    diseaseSpecificFields
    
    # Check everything
    getStatus <- function(x) {
      tmpGenes <- x[[3]]
      tmpGenes <- trimws(strsplit(tmpGenes, ",")[[1]])
      tmpOut <- sapply(tmpGenes, FUN=grepl, x=tmpGeneFindings[,1])
      paste(paste(tmpGeneFindings[as.logical(rowSums(tmpOut)),1], ":", tmpGeneFindings[as.logical(rowSums(tmpOut)),2], sep=""), collapse=", ")
    }
    
    diseaseSpecificFields[,"Value"] <- apply(diseaseSpecificFields, FUN=getStatus, MARGIN=1)
    diseaseSpecificFields[diseaseSpecificFields[,"Value"]==":","Value"] <- "Normal"
    diseaseSpecificFields <- diseaseSpecificFields[,c("Field_name", "Value")]
  } else {
    diseaseSpecificFields <- data.frame()
  }
  return(diseaseSpecificFields)
}


germlineInformation <- function() {
  tmpGeneFindings <- allFindingsTable()
  tmpGeneFindings <- tmpGeneFindings[!grepl("Pathway", tmpGeneFindings[,"Type"]),]
  diseaseSpecificFields

  # Check everything
  getStatus <- function(x) {
    tmpGenes <- x[[3]]
    tmpGenes <- trimws(strsplit(tmpGenes, ",")[[1]])
    tmpOut <- sapply(tmpGenes, FUN=grepl, x=tmpGeneFindings[,1])
    paste(paste(tmpGeneFindings[as.logical(rowSums(tmpOut)),1], ":", tmpGeneFindings[as.logical(rowSums(tmpOut)),2], sep=""), collapse=", ")
  }
  
  diseaseSpecificFields[,"Value"] <- apply(diseaseSpecificFields, FUN=getStatus, MARGIN=1)
  diseaseSpecificFields[diseaseSpecificFields[,"Value"]==":","Value"] <- "Normal"
  diseaseSpecificFields <- diseaseSpecificFields[,c("Field_name", "Value")]
  return(diseaseSpecificFields)
}


genomicSummary <- function() {
  if(exists('expData')){
    tmpRightHead <- c("High Confidence Genomic Alterations", "Total Genomic Alterations", "Transcriptomic Alterations", "Proteomic Alterations", "Aberrant Pathway Activity")
    
    numLesions <- allFindingsTable()
    numLesions <- nrow(numLesions[numLesions[,"Type"]%in%c("Mutation", "Fusion", "Amplification", "Deletion"),])
    numTranscripts <- RNASeqAnalysisOut[[1]][[2]]
    numTranscripts <- nrow(numTranscripts[numTranscripts[,1]>3,])
    numProeins <- "NA"
    numPathways <- RNASeqAnalysisOut[[2]][[2]]
    numPathways <- numPathways[numPathways[,"P_VAL"]<0.01,]
    numPathways <- nrow(numPathways)
    
    tmpVals <- c(nrow(highConfidenceFindingsTable()),numLesions, numTranscripts, numProeins, numPathways)
    df1 <- data.frame(tmpRightHead, tmpVals)
  } else {
    df1 <- data.frame()
  }
  
  return(df1)
}
