# function to convert shared pathways to plots
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyverse)
  library(dplyr)
})

shared_pathways_plot <- function(pathway_analysis, prefix, output_dir){
  shared_pathways <- pathway_analysis$shared_pathways
  
  # if the output for shared pathways is empty
  if(nrow(shared_pathways) == 0){
    return(NULL)
  }
  
  # generate plots for each comparison
  comparisons <- unique(shared_pathways$comparison)
  plist <- list()
  for(i in 1:length(comparisons)){
    shared_pathways_subset <- shared_pathways %>%
      filter(comparison %in% comparisons[i]) %>%
      # arrange by counts
      arrange(desc(sample_count_per_pathway))  %>% 
      select(pathway, sample_count_per_pathway, direction) %>% 
      distinct() %>%
      # take the top 10 for each group
      group_by(direction) %>% 
      dplyr::slice(1:10) %>% 
      arrange(sample_count_per_pathway) %>% 
      mutate(direction = factor(direction, levels = c("up", "down")),
             pathway = factor(pathway, levels = unique(pathway)))
    
    ### Plot for pathway enrichment
    plist[[i]] <- ggplot(shared_pathways_subset, aes(pathway, y = sample_count_per_pathway, fill = direction)) + 
      geom_bar(stat = "identity") + coord_flip() + theme_bw() +
      xlab("") + 
      ylab("Count of Enriched Pathways in 20 Transcriptomically Similar Patients") + 
      theme(plot.margin = unit(c(1, 5, 1, 5), "cm")) + 
      scale_x_discrete(labels = function(x) stringr::str_wrap(gsub('_',' ',x), width = 65)) + 
      scale_fill_manual(name = "Direction", 
                        values = c("up" = "red", "down" = "forest green")) + 
      ggtitle(paste("Comparison vs", comparisons[i]))
  }
  
  fname <- file.path(output_dir, paste0("pathway_analysis_", prefix, ".pdf"))
  pdf(file = fname, width = 12, height = 6)
  for(i in 1:length(plist)){
    print(plist[[i]])
  }
  dev.off()
}
