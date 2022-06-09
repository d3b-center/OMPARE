# barplot of differentially expressed genes

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

diffexpr_genes_barplot <- function(genes_up, genes_down, comparison_study, cancer_genes) {
  
  # filter by comparison
  geneDataUp <- genes_up %>%
    filter(comparison == comparison_study,
           diff_expr == "up",
           genes %in% cancer_genes) %>%
    arrange(desc(logFC)) %>%
    slice_head(n = 20)
  
  geneDataDown <- genes_down %>%
    filter(comparison == comparison_study,
           diff_expr == "down",
           genes %in% cancer_genes) %>%
    arrange(logFC) %>%
    slice_head(n = 20)
  
  # combine both
  geneData <- rbind(geneDataUp, geneDataDown)
  geneData <- geneData %>%
    mutate(Direction = diff_expr,
           Gene = genes) %>%
    arrange(logFC)
  geneData$Gene <- factor(geneData$Gene, levels = geneData$Gene)
  geneData$Direction <- factor(geneData$Direction, levels = c("up", "down"))
  
  # plot barplot of top 20 up/down genes
  p <- ggplot(geneData, aes(Gene, y = logFC, fill = Direction)) + 
    geom_bar(stat="identity") + coord_flip() + theme_bw() + 
    xlab("") + scale_fill_manual(values = c("up" = "red", "down" = "forest green")) +
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
          legend.direction = "none",
          legend.position = "none")  +
    ggtitle(comparison_study) +
    guides(fill = "none")
  return(p)
}
