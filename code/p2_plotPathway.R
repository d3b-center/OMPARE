##########################
# Function to plot pathway
##########################

plotPathway <- function(myRNASeqAnalysisOut = RNASeqAnalysisOut) {
  pathData <- myRNASeqAnalysisOut[[2]][[2]]
  
  # only significant pathways
  pathData <- pathData %>%
    filter(ADJ_P_VALUE < 0.05)
  
  # top 10 upregulated pathways
  pathDataUp <- pathData %>%
    filter(Direction == "Up") %>%
    top_n(10, wt = rev(P_VAL))
  
  # top 10 downregulated pathways
  pathDataDown <- pathData %>%
    filter(Direction == "Down") %>%
    top_n(10, wt = rev(P_VAL))

  # combine and plot
  pathData <- rbind(pathDataDown, pathDataUp)
  pathData$Pathway <- factor(pathData$Pathway, levels = unique(pathData$Pathway))
  pathData$Direction <- factor(pathData$Direction, levels = c("Down", "Up"))
  p <- ggplot(pathData, aes(Pathway, y = (-1)*log10(P_VAL), fill=Direction)) + 
    geom_bar(stat="identity") + coord_flip() + theme_bw() +
    xlab("Pathway Name") + 
    ylab("-log10 P-Value") + scale_fill_manual(values = c("forest green", "red"))
  return(p)
}