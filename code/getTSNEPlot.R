##################################
# Function to generate t-SNE plot
##################################

# Top 10000 most variable genes
getTSNEPlot <- function(dat, clinData, study, patient) {
  set.seed(100) # set seed for reproducibility
  tsneOut <- Rtsne(t(log2(dat+1)), initial_dims=100, perplexity=30, check_duplicates = FALSE)
  tsneData <- data.frame(tsneOut$Y, samples = colnames(dat))
  colnames(tsneData)[1:2] <- c("PC1", "PC2")
  tsneData <- merge(clinData, tsneData, by.x = "sample_barcode", by.y = "samples")
  
  # reverse factors to set PNOC to lowest
  # this is for shape
  tsneData$label <- tsneData$study_id
  tsneData$label <- ifelse(tsneData$sample_barcode == patient, patient, tsneData$label)
  tsneData$label <- factor(tsneData$label, levels = c(study, "PNOC008",patient))

  # this is for size
  tsneData$study_id <- as.factor(tsneData$study_id)
  tsneData$study_id <- relevel(tsneData$study_id, ref = "PNOC008")
  tsneData$study_id <- fct_rev(tsneData$study_id)
  if(study != "PBTA"){
    tsneData$pathology_diagnosis <- tsneData$short_histology
    tsneData$integrated_diagnosis <- tsneData$broad_histology
  }
  
  # plot t-SNE
  p <- ggplot(tsneData, aes(PC1, PC2, 
                            color = pathology_diagnosis,
                            size = study_id,
                            shape = label,
                            text = paste0("Sample:",tsneData$sample_barcode,
                                          "\nShort_histology:", tsneData$short_histology,
                                          "\nBroad_histology:", tsneData$broad_histology,
                                          "\nPathology_diagnosis:", tsneData$pathology_diagnosis,
                                          "\nIntegrated_diagnosis:", tsneData$integrated_diagnosis))) +
    geom_jitter(width = 0.5, height = 0.5) +
    theme_bw() + ggtitle("T-SNE Clustering") +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          legend.title=element_text(size=12), 
          legend.text=element_text(size=12)) + 
    guides(size = FALSE, shape = F, color = F)
  ggplotly(p, tooltip = "text")
}