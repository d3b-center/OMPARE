getTSNEPlot <- function() {
  tsneOut <- Rtsne(t(log2(resTmp+1)), initial_dims=100, perplexity=30, check_duplicates = FALSE)
  tsneData <- data.frame(tsneOut$Y, colnames(resTmp))
  tsneData <- cbind(clinData, tsneData)
  tsneData[,"SampleX"] <- ifelse(tsneData[,"V2"]=="PatSample", 2, 1)
  tmpCol <- as.character(clinData[,"Cancer.Type"])
  p <- ggplot(tsneData, aes(X1, X2, color=Cancer.Type, size=SampleX, shape=as.character(SampleX)))+geom_jitter(width = 0.5, height = 0.5)+theme_bw()+ggtitle("T-SNE PBTA RNA-Sequencing")
  p <- p+theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
               legend.title=element_text(size=12), 
               legend.text=element_text(size=12))+ guides(shape=FALSE, size=FALSE)
  p  
}