# Author: Komal S. Rathi
# Function: Connectivity Analysis on cohort of interest using LINCS

suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(optparse)
})

# source functions
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "mut_distance_calc")
source(file.path(module_dir, "utils", "lincs_connectivity.R"))

option_list <- list(
  make_option(c("--patient"), type = "character",
              help = "sample_id of patient of interest"),
  make_option(c("--nb_cluster_assigned"), type = "character",
              help = "cancer group cluster assignment file from coseq_detect"),
  make_option(c("--vtest_output"), type = "character",
              help = "v-test output file (.RDS)"),
  make_option(c("--wtcs_fdr_cutoff"), type = "character",
              help = "cutoff for WTCS"),
  make_option(c("--trend_val"), type = "character",
              help = "trend value i.e. up or down"),
  make_option(c("--cor_score_cutoff"), type = "character", 
              help = "cutoff for correlation value"),
  make_option(c("--num_sets"), type = "character",
              help = "number of sets for network plots"),
  make_option(c("--output_dir"), type = "character",
              help = "output directory of patient of interest")
)
# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
patient <- opt$patient
nb_cluster_assigned <- opt$nb_cluster_assigned
vtest_output <- opt$vtest_output
wtcs_fdr_cutoff <- as.numeric(opt$wtcs_fdr_cutoff)
trend_val <- opt$trend_val
cor_score_cutoff <- as.numeric(opt$cor_score_cutoff)
num_sets <- as.numeric(opt$num_sets)
output_dir <- opt$output_dir

# define directory
plots_dir <- file.path(output_dir, "plots")
dir.create(path = plots_dir, showWarnings = F, recursive = T)

# find the cluster assigned for the patient
cluster_assigned <- readr::read_tsv(nb_cluster_assigned) %>%
  pull(cluster_assigned_nb)

# filter vtest scores to only cluster of interest 
vtest_output <- readRDS(vtest_output) %>%
  dplyr::filter(cluster == cluster_assigned)

# LINCS-based similarity metric
drug_pathways_barplot <- lincs_connectivity(input = vtest_output,
                                            method = "LINCS",
                                            wtcs_fdr_cutoff = wtcs_fdr_cutoff,
                                            trend_val = trend_val,
                                            cor_score_cutoff = cor_score_cutoff,
                                            output_dir = output_dir)

# save output
fname <- file.path(output_dir, "drug_pathways_barplot.pdf")
pdf(file = fname, width = 18, height = 8)
if(is.null(drug_pathways_barplot)){
  plot.new()
} else {
  drug_pathways_barplot
}
dev.off()
