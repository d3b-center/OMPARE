suppressPackageStartupMessages({
  library(tidyverse)
  library(ggridges)
  library(ggplot2)
  library(dplyr)
})

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# functions
source(file.path(root_dir, "code", "utils", "pubTheme.R"))

# create density plots and save output
drug_dge_density_plots <- function(combined_data, pnoc008_sample, dge_gene){
  vline <- combined_data %>%
    filter(sample == pnoc008_sample) %>%
    dplyr::mutate(tpm = log2(tpm + 1)) %>%
    .$tpm
  combined_data <- combined_data %>%
    filter(!sample %in% c(pnoc008_sample))
  
  # add labels
  combined_data <- combined_data %>%
    group_by(study_id) %>%
    dplyr::mutate(study_id = paste0(study_id, "\n(n = ",n(),")"))
  
  # compute stats
  d2 <- combined_data %>%
    group_by(study_id) %>%
    summarize(lower = quantile(log2(tpm + 1), probs = .025),
              median = quantile(log2(tpm + 1), probs = .5),
              upper = quantile(log2(tpm + 1), probs = .975))
  
  p <- ggplot(data = combined_data, aes(x = log2(tpm + 1))) +
    geom_density(aes(fill = study_id, alpha = 0.5)) +
    theme_Publication(base_size = 10) +
    facet_wrap(~study_id, scales = "free") +
    geom_vline(xintercept = vline, color = "red") + 
    geom_vline(data = d2, aes(xintercept = upper), color = "black", linetype = "dashed") + 
    geom_vline(data = d2, aes(xintercept = median), color = "black", linetype = "dashed") + 
    geom_vline(data = d2, aes(xintercept = lower), color = "black", linetype = "dashed") + 
    guides(alpha = F,  fill = F) +
    ggtitle(dge_gene) +
    xlab("log2 TPM") 
  
  return(p)
}

