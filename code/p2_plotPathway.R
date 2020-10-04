##########################
# Function to plot pathway
##########################

plotPathway <- function(myRNASeqAnalysisOut = RNASeqAnalysisOut) {
  pathData <- myRNASeqAnalysisOut$pathways
  
  # only significant pathways
  pathData <- pathData %>%
    filter(ADJ_P_VAL < 0.05)
  
  # top 10 upregulated pathways
  pathDataUp <- pathData %>%
    filter(Direction == "Up") %>%
    arrange(ADJ_P_VAL) %>%
    slice_head(n = 10)
  
  # top 10 downregulated pathways
  pathDataDown <- pathData %>%
    filter(Direction == "Down") %>%
    arrange(ADJ_P_VAL) %>%
    slice_head(n = 10)

  # combine and plot
  pathData <- rbind(pathDataDown, pathDataUp)
  pathData$Pathway <- factor(pathData$Pathway, levels = unique(pathData$Pathway))
  pathData$Direction <- factor(pathData$Direction, levels = c("Down", "Up"))
  p <- ggplot(pathData, aes(Pathway, y = (-1)*log10(ADJ_P_VAL), fill=Direction)) + 
    geom_bar(stat="identity") + coord_flip() + theme_bw() +
    xlab("") + 
    ylab("-log10 Adj. P-Value") + scale_fill_manual(values = c("forest green", "red")) +
    theme(plot.margin = unit(c(1, 5, 1, 7), "cm"))
  return(p)
}