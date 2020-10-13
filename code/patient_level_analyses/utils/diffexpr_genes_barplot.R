# barplot of differentially expressed genes

diffexpr_genes_barplot <- function(rnaseq_analysis_output = rnaseq_analysis_output) {
  geneData <- rnaseq_analysis_output$diffexpr.top20 %>%
    rownames_to_column("Gene") %>%
    mutate(Direction = ifelse(z_score > 0, "Up", "Down")) %>%
    arrange(z_score)
  geneData$Gene <- factor(geneData$Gene, levels = geneData$Gene)
  
  # plot barplot of top 20 up/down genes
  p <- ggplot(geneData, aes(Gene, y = z_score, fill = Direction)) + 
    geom_bar(stat="identity") + coord_flip() + theme_bw() + 
    xlab("") + scale_fill_manual(values = c("Down" = "forest green", "Up" = "red")) +
    theme(plot.margin = unit(c(1, 5, 1, 5), "cm"))
  return(p)
}