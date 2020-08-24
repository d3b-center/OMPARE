# Genomic Summary
genomicSummary <- function() {
  if(exists('expData')){
    headers <- c("High Confidence Genomic Alterations", "Total Genomic Alterations", "Transcriptomic Alterations", "Aberrant Pathway Activity")
    
    # high confidence alterations
    highConfLesions <- nrow(highConfidenceFindingsTable())
    numLesions <- allFindingsTable() %>%
      filter(Type %in% c("Mutation", "Fusion", "Amplification", "Deletion")) %>%
      nrow()
    
    # highly upreg genes (z-score > 3)
    numTranscripts <- RNASeqAnalysisOut[[1]][[2]] %>%
      filter(Z_Score > 3) %>%
      nrow()
    
    # adj. pval < 0.05 (highly significant pathways)
    numPathways <- RNASeqAnalysisOut[[2]][[2]] %>%
      filter(ADJ_P_VALUE < 0.05) %>%
      nrow()
    
    tmpVals <- c(highConfLesions, numLesions, numTranscripts, numPathways)
    df1 <- data.frame(headers, tmpVals)
  } else {
    df1 <- data.frame()
  }
  
  return(df1)
}