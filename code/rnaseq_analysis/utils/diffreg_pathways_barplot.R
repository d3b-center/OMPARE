# barplot of differentially regulated up/down pathways
diffreg_pathways_barplot <- function(pathways_up, pathways_down, comparison_study) {
  
  # top 10 upregulated pathways
  pathDataUp <- pathways_up %>%
    filter(comparison == comparison_study) %>%
    arrange(padj) %>%
    slice_head(n = 10)
  
  # top 10 downregulated pathways
  pathDataDown <- pathways_down %>%
    filter(comparison == comparison_study) %>%
    arrange(padj) %>%
    slice_head(n = 10)
  
  # combine and plot
  pathData <- rbind(pathDataDown, pathDataUp)
  pathData$direction <- str_to_title(pathData$direction)
  pathData$pathway <- factor(pathData$pathway, levels = unique(pathData$pathway))
  pathData$direction <- factor(pathData$direction, levels = c("Up", "Down"))
  p <- ggplot(pathData, aes(pathway, y = (-1)*log10(padj), fill = direction)) + 
    geom_bar(stat="identity") + coord_flip() + theme_bw() +
    xlab("") + 
    ylab("-log10 Adj. P-Value") + 
    scale_fill_manual(name = "Direction", values = c("Down" = "forest green", "Up" = "red")) +
    theme(plot.margin = unit(c(1, 5, 1, 7), "cm")) + 
    ggtitle(paste0("Comparison against ", comparison_study))
  return(p)
}
