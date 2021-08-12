library(tidyverse)
library(dplyr)
library(ggplot2)

plot_immune_profile  <- function(xcell_scores){
  # format data
  xcell_scores <- xcell_scores %>%
    rownames_to_column('CellType') %>%
    gather("Sample", "Score", -CellType) %>%
    mutate("IsSample" = ifelse(grepl(sampleInfo$subjectID, Sample), T, F))
  xcell_scores.sample <- xcell_scores %>%
    filter(IsSample == TRUE)
  
  # set factors
  celltype.order <- xcell_scores %>%
    group_by(CellType) %>%
    summarise(median = median(Score)) %>%
    arrange(desc(median)) %>%
    .$CellType
  xcell_scores$CellType <- factor(xcell_scores$CellType, levels = celltype.order)
  
  # boxplot
  p <- ggplot(xcell_scores, aes(CellType, Score)) + 
    geom_boxplot(outlier.shape = NA) +  
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
    geom_point(data = xcell_scores.sample, aes(CellType, Score), colour = "red", size = 2, shape = "triangle") +
    theme(axis.text = element_text(size = 8, face = "bold"), 
          axis.title = element_blank())
  return(p)
}
