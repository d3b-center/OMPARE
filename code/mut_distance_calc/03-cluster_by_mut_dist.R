# Author: Run Jin
# Obtain NB.MClust assignment of POI based on 20 closest neighbors

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(readr)
})

#### Parse command line options ------------------------------------------------
option_list <- list(
  make_option(c("--ref_cancer_dir"),type="character",
              help="directory with cancer type specific reference tumor subset"),
  make_option(c("--sample_id_interest"), type = "character",
              help = "sample_id of patient of interest"),
  make_option(c("--coseq_cluster_assignment"), type = "character",
              help = "path to coseq cluster assignment file"),
  make_option(c("--mut_dist"), type = "character",
              help = "file containing mutational distance between poi vs. cohort (.tsv)"),
  make_option(c("--output_dir"), type = "character",
              help = "output directory of patient of interest")
)
opt <- parse_args(OptionParser(option_list=option_list))
ref_cancer_dir <- opt$ref_cancer_dir
sample_id_interest <- opt$sample_id_interest
coseq_cluster_assignment <- opt$coseq_cluster_assignment
mut_dist <- opt$mut_dist
output_dir <- opt$output_dir

# Read in files necessary for analyses
histology_df <- list.files(path = ref_cancer_dir, pattern = "histologies", full.names = T)
histology_df <- readr::read_tsv(histology_df)

# cluster assignment file
coseq_cluster_assignment <- readr::read_tsv(coseq_cluster_assignment)

# mutational distance file
mut_dist <- readr::read_tsv(mut_dist)

# Find the sample id of the samples in NB.MClust file
sample_id_df <- histology_df %>%
  dplyr::filter(Kids_First_Biospecimen_ID %in% coseq_cluster_assignment$Kids_First_Biospecimen_ID) %>%
  select(sample_id, Kids_First_Biospecimen_ID)

# add the sample ID to the results for mapping 
coseq_cluster_assignment <- coseq_cluster_assignment %>%
  left_join(sample_id_df)

# find the closest 20 samples from mutational distance matrix
closest_20 <- mut_dist %>% 
  arrange(all_of(gower_dist_poi)) %>% 
  slice_head(n = 20) %>%
  pull(sample_id) 

# save closest 20 nearest neighbors + additional info for downstream analysis
nn_table <- mut_dist %>% 
  left_join(histology_df %>%
              dplyr::select(sample_id, OS_days, OS_status)) %>%
  unique() %>%
  dplyr::mutate(nearest_neighbor = sample_id,
                distance = gower_dist_poi)
nn_table %>% 
  filter(sample_id %in% closest_20) %>%
  dplyr::select(nearest_neighbor, distance) %>%
  saveRDS(file.path(output_dir, "nn_table.rds"))
nn_table %>% 
  dplyr::select(sample_id, OS_days, OS_status) %>%
  saveRDS(file.path(output_dir, "surv_data.rds"))

# get the cluster of the 20 closest samples 
count_assignment <- coseq_cluster_assignment %>% 
  filter(sample_id %in% closest_20) %>% 
  group_by(cluster_assigned_nb ) %>%
  dplyr::mutate(n=n()) 

max_count <- max(count_assignment$n)
cluster_assigned <- count_assignment %>% 
  dplyr::filter(n == max_count) %>% 
  dplyr::select(cluster_assigned_nb) %>%
  distinct()

################ Finally output the cluster number of sample of interest
cluster_assigned %>% readr::write_tsv(file.path(output_dir, "patient_of_interest_nb_cluster_assigned.tsv"))

