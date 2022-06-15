# Author: Run Jin
# Calculate distance of samples in the cohort by measurements of mutations

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(readr)
  library(StatMatch)
})

#### Parse command line options ------------------------------------------------
option_list <- list(
  make_option(c("--sample_id_interest"), type = "character",
              help = "`sample_id` of patient of interest"),
  make_option(c("--mat_file"), type = "character",
              help = "binary matrix file with gene vs. sample id indicating mutational status (.tsv)"),
  make_option(c("--dist_file"), type = "character",
              help = "file containing the mutational distance results (.tsv)")
)
opt <- parse_args(OptionParser(option_list=option_list))
sample_id_interest <- opt$sample_id_interest
mat_file <- opt$mat_file
dist_file <- opt$dist_file

# read in files necessary for analyses
# this file already has POI so no need to calculate POI binary matrix
mut_matrix <- readr::read_tsv(mat_file) %>%
  tibble::column_to_rownames("...1")

# calculate distance
dx <- gower.dist(t(as.matrix(mut_matrix)))
dx <- dx[sample_id_interest,] %>% as.data.frame()

# extract distance
colnames(dx) <- "gower_dist_poi"

# write output
dx %>% tibble::rownames_to_column("sample_id") %>%
  dplyr::filter(sample_id != sample_id_interest) %>% 
  readr::write_tsv(dist_file)
