# barplot of differentially expressed genes

diffexpr_genes_barplot <- function(rnaseq_analysis_output = rnaseq_analysis_output) {
  
  geneData <- rnaseq_analysis_output$diffexpr.top20 %>%
    rownames_to_column("Gene") %>%
    mutate(Direction = ifelse(logfc > 0, "Up", "Down")) %>%
    arrange(logfc)
  geneData$Gene <- factor(geneData$Gene, levels = geneData$Gene)
  
  # plot barplot of top 20 up/down genes
  geneData$Direction <- factor(geneData$Direction, levels = c("Up", "Down"))
  p <- ggplot(geneData, aes(Gene, y = logfc, fill = Direction)) + 
    geom_bar(stat="identity") + coord_flip() + theme_bw() + 
    xlab("") + scale_fill_manual(values = c("Down" = "forest green", "Up" = "red")) +
    theme(plot.margin = unit(c(1, 5, 1, 5), "cm"))
  return(p)
}