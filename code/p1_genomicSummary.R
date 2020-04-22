# Genomic Summary
genomicSummary <- function() {
  if(exists('expData')){
    tmpRightHead <- c("High Confidence Genomic Alterations", "Total Genomic Alterations", "Transcriptomic Alterations", "Proteomic Alterations", "Aberrant Pathway Activity")
    
    highConfLesions <- nrow(highConfidenceFindingsTable())
    numLesions <- allFindingsTable()
    numLesions <- nrow(numLesions[numLesions[,"Type"]%in%c("Mutation", "Fusion", "Amplification", "Deletion"),])
    numTranscripts <- RNASeqAnalysisOut[[1]][[2]]
    numTranscripts <- nrow(numTranscripts[numTranscripts$Z_Score>3,]) # z-score > 3 (highly upreg genes)
    numProteins <- "NA"
    numPathways <- RNASeqAnalysisOut[[2]][[2]]
    numPathways <- numPathways[numPathways[,"P_VAL"]<0.01,]  # pval < 0.01 (highly significant pathways)
    numPathways <- nrow(numPathways)
    
    tmpVals <- c(highConfLesions,numLesions, numTranscripts, numProteins, numPathways)
    df1 <- data.frame(tmpRightHead, tmpVals)
  } else {
    df1 <- data.frame()
  }
  
  return(df1)
}