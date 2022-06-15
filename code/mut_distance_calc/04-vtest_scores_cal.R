# Author: Run Jin
# Obtain V-test scores for cluster of interest

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(readr)
  library(edgeR)
})

# source functions
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "mut_distance_calc")
source(file.path(module_dir, "utils", "vtest_cal.R"))

#### Parse command line options ------------------------------------------------
option_list <- list(
  make_option(c("--ref_cancer_dir"),type="character",
              help="directory with cancer type specific reference tumor subset"),
  make_option(c("--coseq_cluster_assignment"), type = "character",
              help = "path to coseq cluster assignment file"),
  make_option(c("--filtered_count_cg_coding_pval"), type = "character",
              help = "path to filtered protein coding counts file"),
  make_option(c("--vtest_file"), type = "character",
              help = "output file of vtest scores (.rds)")
)

opt <- parse_args(OptionParser(option_list=option_list))
ref_cancer_dir <- opt$ref_cancer_dir
coseq_cluster_assignment <- opt$coseq_cluster_assignment
filtered_count_cg_coding_pval <- opt$filtered_count_cg_coding_pval
output_dir <- opt$output_dir
vtest_file <- opt$vtest_file

# Read in files necessary for analyses
histology_df <- list.files(path = ref_cancer_dir, pattern = "histologies", full.names = T)
histology_df <- readr::read_tsv(histology_df)

exp_count_cg_coding <- list.files(path = ref_cancer_dir, pattern = "count", full.names = T)
exp_count_cg_coding <- readRDS(exp_count_cg_coding)

# cluster assignment file
coseq_cluster_assignment <- readr::read_tsv(coseq_cluster_assignment)

# read in filtered count cg coding pval file from `coseq_detect.R` output 
filtered_count_cg_coding_pval <- readRDS(filtered_count_cg_coding_pval)

filtered_count_cg_coding_pval_cluster <- filtered_count_cg_coding_pval %>% t() %>% 
  as.data.frame() %>%
  tibble::rownames_to_column('sample') %>%
  tidyr::gather(gene_symbol, gene_count, -sample) %>%
  dplyr::left_join(coseq_cluster_assignment %>%
                     dplyr::rename(sample = Kids_First_Biospecimen_ID) %>%
                     select(sample, cluster_assigned_nb), by = 'sample') 

# gather information about 
filtered_count_cg_coding_pval_cluster <- filtered_count_cg_coding_pval_cluster %>%
  dplyr::group_by(cluster_assigned_nb, gene_symbol) %>%
  dplyr::mutate(cluster_gene_mean_score = mean(gene_count)) %>% # mean of gene count per cluster 
  ungroup() %>%
  dplyr::group_by(gene_symbol) %>%
  dplyr::mutate(gene_mean_score = mean(gene_count),
                gene_variance = var(gene_count)) # global mean & variance per gene count

# apply v.test function per gene
vtest_output <- plyr::ddply(.data = filtered_count_cg_coding_pval_cluster, 
                            .variables = "gene_symbol", 
                            .fun = function(x) compute.v.test(x, clustering_col = "cluster_assigned_nb"))

vtest_output <- vtest_output %>%
  dplyr::rename("geneSymbol" = "gene_symbol", "score" = "v_score") %>% 
  saveRDS(vtest_file)


