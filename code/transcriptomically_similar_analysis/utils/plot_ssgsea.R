suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(ggplot2)
})

plot_ssgsea <- function(ssgsea_input) { 
  # factorize by median
  tmp <- ssgsea_input %>%
    dplyr::select(geneset_name, gsea_score_median) %>%
    arrange(desc(gsea_score_median)) %>%
    unique() %>%
    .$geneset_name
  ssgsea_input$geneset_name <- factor(ssgsea_input$geneset_name, levels = tmp)
  raw.scoresSample <- ssgsea_input[ssgsea_input$IsSample == T,]
  
  # plot as boxplot
  p <- ggplot(ssgsea_input, aes(geneset_name, gsea_score)) + 
    geom_boxplot() +  
    theme_bw(base_size = 10) +
    geom_point(data = raw.scoresSample, aes(geneset_name, gsea_score), colour = "red") +
    theme(axis.title = element_blank()) + coord_flip()

  return(p)  
}
