# barplot of differentially expressed genes

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

diffexpr_genes_barplot <- function(genes_diffexpr, cancer_genes) {
  
  # filter by comparison
  geneDataUp <- genes_diffexpr %>%
    filter(diff_expr == "up",
           genes %in% cancer_genes) %>%
    arrange(desc(logFC)) %>%
    slice_head(n = 20)
  
  geneDataDown <- genes_diffexpr %>%
    filter(diff_expr == "down",
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
  
  # comparison study
  comparison_study <- unique(geneData$comparison)
  
  # plot barplot of top 20 up/down genes
  p <- ggplot(geneData, aes(Gene, y = logFC, fill = Direction)) + 
    geom_bar(stat="identity") + coord_flip() + theme_bw() + 
    xlab("") + scale_fill_manual(values = c("up" = "red", "down" = "forest green")) +
    theme(legend.direction = "none",
          legend.position = "none")  +
    ggtitle(paste0("Comparison: ", comparison_study)) +
    guides(fill = "none")
  return(p)
}
