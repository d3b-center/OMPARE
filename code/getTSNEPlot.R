getTSNEPlot <- function() {
  set.seed(100) # set seed for reproducibility
  tsneOut <- Rtsne(t(log2(resTmp+1)), initial_dims=100, perplexity=30, check_duplicates = FALSE)
  tsneData <- data.frame(tsneOut$Y, samples = colnames(resTmp))
  colnames(tsneData)[1:2] <- c("PC1", "PC2")
  tsneData <- merge(clinData, tsneData, by.x = "Kids_First_Biospecimen_ID", by.y = "samples")
  tsneData$type <- ifelse(tsneData$sample_id==sampleInfo$subjectID, "PNOC", "PBTA")
  tsneData$type <- factor(tsneData$type)
  
  p <- ggplot(tsneData, aes(PC1, PC2, 
                       color = pathology_diagnosis,
                       size = type,
                       shape = type,
                       text = paste0("Participant:",tsneData$Kids_First_Biospecimen_ID,
                                     "\nSample:", tsneData$sample_id,
                                     "\nShort_histology:", tsneData$short_histology,
                                     "\nBroad_histology:", tsneData$broad_histology,
                                     "\nPathology_diagnosis:", tsneData$pathology_diagnosis,
                                     "\nIntegrated_diagnosis:", tsneData$integrated_diagnosis))) +
    geom_jitter(width = 0.5, height = 0.5) +
    theme_bw() + ggtitle("T-SNE PBTA RNA-Sequencing") +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          legend.title=element_text(size=12), 
          legend.text=element_text(size=12)) + 
    guides(size = FALSE, color = F)
  ggplotly(p, tooltip = "text")
}
