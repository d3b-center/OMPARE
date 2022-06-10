suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
})

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# output directory
module_dir <- file.path(root_dir, "code", "drug_recommendations")
output_dir <- file.path(patient_dir, "output", "drug_recommendations")
output_dir <- file.path(output_dir, "drug_dge_density_plots")
dir.create(output_dir, showWarnings = F, recursive = T)

# source functions
source(file.path(module_dir, "utils", "drug_dge_density_plots.R"))

# input data
transcriptomic_drug_rec <- readRDS(file.path(patient_dir, "output", "drug_recommendations", "transcriptome_drug_rec.rds"))
dge_all <- transcriptomic_drug_rec %>% 
  dplyr::select(Gene, Comparison, logFC) %>%
  unique()

patient_of_interest <- sample_info %>%
  filter(experimental_strategy == "RNA-Seq") %>%
  pull(Kids_First_Biospecimen_ID)

# call function to generate drug density plots
drug_dge_density_plots(patient_of_interest = patient_of_interest, 
                       dge_all = dge_all,
                       output_dir = output_dir)

