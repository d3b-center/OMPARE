library(tidyverse)
library(dplyr)
library(ggplot2)

plot_ssgsea <- function(ssgsea_input) { 
  # factorize by median
  tmp <- ssgsea_input %>%
    dplyr::select(geneset_name, gsea_score_median) %>%
    arrange(desc(gsea_score_median)) %>%
    unique() %>%
    .$geneset_name
  ssgsea_input$geneset_name <- factor(ssgsea_input$geneset_name, levels = tmp)
  
  # plot as boxplot
  p <- ggplot(ssgsea_input, aes(geneset_name, gsea_score)) + 
    geom_boxplot(outlier.shape = NA) +  
    theme_bw()
  raw.scoresSample <- ssgsea_input[ssgsea_input$IsSample == T,]
  p <- p + 
    geom_point(data = raw.scoresSample, aes(geneset_name, gsea_score), colour = "red", size = 3, shape = "triangle") +
    theme(axis.text = element_text(size = 8, face = "bold"), 
          axis.title = element_blank()) + coord_flip() +
    theme(plot.margin = unit(c(1, 5, 1, 2), "cm"))
  return(p)  
}
