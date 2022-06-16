# function to generate fusion data 
suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
})

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
module_dir <- file.path(root_dir, "code", "create_background_matrices")

# source functions
source(file.path(module_dir, "utils", "filter_fusions.R"))

# cancer genes 
cancer_genes <- readRDS(file.path(data_dir, "cancer_gene_list.rds"))

create_fusion_files <- function(star_fusion, arriba_fusion, hist_file, output_dir, prefix){
  
  # create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # subset to id columns
  hist_file <- hist_file %>%
    dplyr::select(Kids_First_Biospecimen_ID, sample_id, cohort, cohort_participant_id) %>%
    unique()
  
  # filter fusions
  star_fusion <- filter_fusions(fusion_data = star_fusion, myCancerGenes = cancer_genes, method = "star-fusion")
  star_fusion$method <- "star-fusion"
  arriba_fusion <- filter_fusions(fusion_data = arriba_fusion, myCancerGenes = cancer_genes, method = "arriba")
  arriba_fusion$method <- "arriba_fusion"
  fusion_data <- rbind(star_fusion, arriba_fusion)
  
  # join with histology
  fusion_data <- fusion_data %>%
    dplyr::rename("Kids_First_Biospecimen_ID" = "tumor_id") %>%
    inner_join(hist_file, by = "Kids_First_Biospecimen_ID") %>%
    group_by(fusion_name, gene1, gene2, annots, Kids_First_Biospecimen_ID, sample_id, cohort, cohort_participant_id) %>%
    summarise(method = toString(method),
              type = toString(type)) %>%
    ungroup()
  
  saveRDS(fusion_data, file = file.path(output_dir, paste(prefix, 'fusion_filtered.rds', sep = "_")))
}
