# barplot of differentially expressed genes

diffexpr_genes_barplot <- function(rnaseq_analysis_output = rnaseq_analysis_output) {
  geneData <- rnaseq_analysis_output$diffexpr.top20
  geneData <- geneData %>%
    rownames_to_column("Gene") %>%
    mutate(Direction = ifelse(Z_Score > 0, "Up", "Down")) %>%
    arrange(Z_Score)
  geneData$Gene <- factor(geneData$Gene, levels = geneData$Gene)
  
  # plot barplot of top 20 up/down genes
  p <- ggplot(geneData, aes(Gene, y = Z_Score, fill = Direction)) + 
    geom_bar(stat="identity") + coord_flip() + theme_bw() + 
    xlab("") + scale_fill_manual(values = c("Down" = "forest green", "Up" = "red")) +
    theme(plot.margin = unit(c(1, 5, 1, 5), "cm"))
  return(p)
}