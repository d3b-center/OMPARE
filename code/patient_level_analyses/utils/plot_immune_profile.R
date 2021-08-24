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
    geom_boxplot() +
    theme_bw(base_size = 10) + 
    geom_point(data = xcell_scores.sample, aes(CellType, Score), colour = "red") +
    theme(axis.title = element_blank()) + coord_flip()
  return(p)
}
