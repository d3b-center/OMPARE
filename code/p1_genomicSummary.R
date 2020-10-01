# Genomic Summary
genomicSummary <- function(snv_pattern) {
  if(exists('expData')){
    headers <- c("High Confidence Genomic Alterations", "Total Genomic Alterations", "Transcriptomic Alterations", "Aberrant Pathway Activity")
    
    # high confidence alterations
    highConfLesions <- nrow(highConfidenceFindingsTable(snv_pattern))
    numLesions <- allFindingsTable(snv_pattern) %>%
      filter(Type %in% c("Mutation", "Fusion", "Amplification", "Deletion")) %>%
      nrow()
    
    # highly upreg genes (z-score > 3)
    numTranscripts <- RNASeqAnalysisOut$diffexpr.top20 %>%
      filter(Z_Score > 3) %>%
      nrow()
    
    # adj. pval < 0.05 (highly significant pathways)
    numPathways <- RNASeqAnalysisOut$pathways %>%
      filter(ADJ_P_VAL < 0.05) %>%
      nrow()
    
    tmpVals <- c(highConfLesions, numLesions, numTranscripts, numPathways)
    df1 <- data.frame(headers, tmpVals)
  } else {
    df1 <- data.frame()
  }
  
  return(df1)
}