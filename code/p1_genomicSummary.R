# Genomic Summary
genomicSummary <- function() {
  if(exists('expData')){
    headers <- c("High Confidence Genomic Alterations", "Total Genomic Alterations", "Transcriptomic Alterations", "Aberrant Pathway Activity")
    
    highConfLesions <- nrow(highConfidenceFindingsTable())
    numLesions <- allFindingsTable()
    numLesions <- nrow(numLesions[numLesions[,"Type"]%in%c("Mutation", "Fusion", "Amplification", "Deletion"),])
    numTranscripts <- RNASeqAnalysisOut[[1]][[2]]
    numTranscripts <- nrow(numTranscripts[numTranscripts$Z_Score>3,]) # z-score > 3 (highly upreg genes)
    numPathways <- RNASeqAnalysisOut[[2]][[2]]
    numPathways <- numPathways[numPathways$ADJ_P_VAL < 0.05,]  # adj. pval < 0.05 (highly significant pathways)
    numPathways <- nrow(numPathways)
    
    tmpVals <- c(highConfLesions, numLesions, numTranscripts, numPathways)
    df1 <- data.frame(headers, tmpVals)
  } else {
    df1 <- data.frame()
  }
  
  return(df1)
}