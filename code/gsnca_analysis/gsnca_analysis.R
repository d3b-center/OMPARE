# Author: Run Jin
# GSNCA analysis comparing upper and lower quantile of gene expressions in each disease

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
})

# arguments
option_list <- list(
  make_option(c("--patient"), type = "character",
              help = "Kids_First_Biospecimen_ID of patient of interest"),
  make_option(c("--pediatric_cancer_dir"),type="character",
              help="directory with cancer type specific pediatric tumor subset"),
  make_option(c("--adult_data_dir"),type="character",
              help="directory to adult data matching specific pediatric tumor subset"),
  make_option(c("--norm_data_dir"), type = "character",
              help = "Cancer group of the patient of interest"),
  make_option(c("--nb_cluster_assigned"), type = "character",
              help = "cancer group cluster assignment file from coseq_detect"),
  make_option(c("--output_dir"), type = "character",
              help = "output directory of patient of interest")
  
)
opt <- parse_args(OptionParser(option_list = option_list))
patient <- opt$patient
pediatric_cancer_dir <- opt$pediatric_cancer_dir
norm_data_dir <- opt$norm_data_dir
adult_data_dir <- opt$adult_data_dir
nb_cluster_assigned <- readr::read_tsv(opt$nb_cluster_assigned)
output_dir <- opt$output_dir

#### Define Directories --------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "gsnca_analysis")

#### Source utils functions --------------------------------------------------------
source(file.path(module_dir, "utils", "gsnca_calc.R"))

#### Read in files necessary for analyses --------------------------------------
# NB.MClust results
cluster_assigned <- nb_cluster_assigned %>%
  filter(Kids_First_Biospecimen_ID == patient) %>%
  pull(cluster_assigned_nb)

# get Kids First Biospecimen IDs of those samples assigned
samples_in_cluster_assigned <- nb_cluster_assigned %>%
  dplyr::filter(cluster_assigned_nb == cluster_assigned) %>%
  pull(Kids_First_Biospecimen_ID)

# Normal Sample TPM
norm_data_tpm <- list.files(path = norm_data_dir, pattern = "tpm", full.names = T)
norm_data_tpm  <- readRDS(norm_data_tpm)

# expression (TPM) of cancer group of interest
cancer_group_tpm <- list.files(path = pediatric_cancer_dir, pattern = "tpm", full.names = T)
cancer_group_tpm <- readRDS(cancer_group_tpm)

# expression (TPM) of cancer group of interest in adult cohort 
adult_data_tpm <- list.files(path = adult_data_dir, pattern = "tpm", full.names = T)
adult_data_tpm <- readRDS(adult_data_tpm)

# expression TPM matrix of MB.MClust samples of cluster of interest 
stopifnot(all_of(samples_in_cluster_assigned) %in% colnames(cancer_group_tpm))
cluster_samples_tpm <- cancer_group_tpm %>% 
  dplyr::select(all_of(samples_in_cluster_assigned)) %>%
  filter_low_expr_df()

#### Prepare background/cohort comparison expression matrix -------------------------------
# normal matrix filtered
norm_data_tpm_filtered <- filter_low_expr_df(norm_data_tpm)

# the rest of the cancer group of interest 
cancer_group_tpm_rest_filtered <- cancer_group_tpm %>%
  dplyr::select(-all_of(samples_in_cluster_assigned)) %>%
  filter_low_expr_df()

# adult matrix filtered 
adult_data_tpm_filtered <- filter_low_expr_df(adult_data_tpm)

#### Run GSNCA and output results in text files ----------------------
gsnca_analysis_plot(cluster_samples_tpm_df = cluster_samples_tpm, 
                    ref_expr_df = norm_data_tpm_filtered, 
                    ref_name = "normal_data", 
                    top_bar=20, 
                    top_net=5,
                    output_dir = output_dir)

gsnca_analysis_plot(cluster_samples_tpm_df = cluster_samples_tpm, 
                    ref_expr_df = cancer_group_tpm_rest_filtered, 
                    ref_name = "pediatric_data", 
                    top_bar=20, 
                    top_net=5,
                    output_dir = output_dir)

gsnca_analysis_plot(cluster_samples_tpm_df = cluster_samples_tpm, 
                    ref_expr_df = adult_data_tpm_filtered, 
                    ref_name = "adult_data", 
                    top_bar=20, 
                    top_net=5,
                    output_dir = output_dir)


