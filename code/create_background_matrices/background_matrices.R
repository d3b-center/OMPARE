# generate rnaseq and non-rnaseq matrices for background comparators 
# adult tumors, pediatric tumors and normal tissues

# load libraries
library(optparse)

# parse params
option_list <- list(
  make_option(c("--normal_tissue"), type = "character",
              help = "normal tissue type"),
  make_option(c("--pediatric_cancer"), type = "character",
              help = "pediatric cancer type"),
  make_option(c("--adult_cancer"), type = "character",
              help = "adult cancer type")
)
opt <- parse_args(OptionParser(option_list = option_list))
normal_tissue <- opt$normal_tissue
pediatric_cancer <- opt$pediatric_cancer
adult_cancer <- opt$adult_cancer

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir('.git'))
data_dir <- file.path(root_dir, "data")
module_dir <- file.path(root_dir, "code", "create_background_matrices")

# source function
source(file.path(module_dir, "utils", "create_histologies_subset.R"))
source(file.path(module_dir, "utils", "create_rnaseq_matrices.R"))
source(file.path(module_dir, "utils", "create_non_rnaseq_matrices.R"))
source(file.path(module_dir, "utils", "create_tcga_non_rnaseq_matrices.R"))
source(file.path(module_dir, "utils", "collapse_rnaseq.R"))
source(file.path(module_dir, "utils", "create_fusion_files.R"))
source(file.path(module_dir, "utils", "diff_expr.R"))

# output directories
normal_tissue_dir <- file.path(data_dir, "normal_data")
pediatric_cancer_dir <- file.path(data_dir, "pediatric_data")
adult_cancer_dir <- file.path(data_dir, "adult_data")

# minimum set of columns used from histology file
# this will help minimize hardcoding for patients not being added to the histology
hist_cols_used <- c("Kids_First_Biospecimen_ID","sample_id","cohort","cohort_participant_id",
                    "aliquot_id", "experimental_strategy","RNA_library",
                    "sample_type","tumor_descriptor",
                    "pathology_diagnosis","integrated_diagnosis",
                    "short_histology","broad_histology","molecular_subtype",
                    "gtex_group",
                    "OS_days","OS_status")

# load master genomics (this will have +1 data)
master_genomics_dir <- file.path(data_dir, "master_genomics")

# histology file
hist_file <- readr::read_tsv(file.path(master_genomics_dir, "histologies-base.tsv"))
hist_file <- hist_file %>%
  dplyr::select(hist_cols_used)

# OT for normal tissues and adult cancers
open_pedcan_analysis_data <- file.path(data_dir, "OpenPedCan-analysis", "data")

# tpm matrices
tpm_mat <- file.path(open_pedcan_analysis_data, "gene-expression-rsem-tpm-collapsed.rds")
tpm_mat <- readRDS(tpm_mat)

# count matrices
counts_mat <- file.path(open_pedcan_analysis_data, "gene-counts-rsem-expected_count-collapsed.rds")
counts_mat <- readRDS(counts_mat)

# 1. normal tissue (from OpenPedCan)

# create histology subset
gtex_hist_file <- create_histologies_subset(hist_file = hist_file, 
                                            group_filter = normal_tissue,
                                            cohort_filter = "GTEx",
                                            output_dir = normal_tissue_dir,
                                            prefix = "normal_tissues")

# normal rnaseq matrices
create_rnaseq_matrices(hist_file = gtex_hist_file, 
                       tpm_mat = tpm_mat, 
                       counts_mat = counts_mat, 
                       output_dir = normal_tissue_dir, 
                       prefix = "normal_tissues")

# 2. adult cancer comparator (from OpenPedCan)

# create histology subset
adult_hist_file <- create_histologies_subset(hist_file = hist_file,
                                             group_filter = adult_cancer,
                                             cohort_filter = "TCGA",
                                             output_dir = adult_cancer_dir,
                                             prefix = "adult_tumors")

# tpm matrices
tpm_mat <- file.path(open_pedcan_analysis_data, "tcga-gene-expression-rsem-tpm-collapsed.rds")
tpm_mat <- readRDS(tpm_mat)

# count matrices
counts_mat <- file.path(open_pedcan_analysis_data, "tcga-gene-counts-rsem-expected_count-collapsed.rds")
counts_mat <- readRDS(counts_mat)

# adult rnaseq matrices
create_rnaseq_matrices(hist_file = adult_hist_file,
                       tpm_mat = tpm_mat,
                       counts_mat = counts_mat,
                       output_dir = adult_cancer_dir,
                       prefix = "adult_tumors")

# adult non-rnaseq matrices (use TCGAbiolinks)
create_tcga_non_rnaseq_matrices(hist_file = adult_hist_file,
                                cohort_filter = "TCGA",
                                group_filter = adult_group_filter,
                                output_dir = adult_cancer_dir,
                                prefix = "adult_tumors")

# 3. pediatric cancer comparator (from master genomics)
# master genomics files
# master_genomics_dir <- file.path(data_dir, "master_genomics")
master_genomics_files <- list.files(path = master_genomics_dir, full.names = T)

# histology file from data assembly project  
master_genomics_hist_file <- master_genomics_files[grep("histologies-base.tsv", master_genomics_files)]  %>%
  readr::read_tsv()

# subset to pediatric cohorts
master_genomics_hist_file <- master_genomics_hist_file %>%
  filter(!cohort %in% c("TCGA", "GTEx"))
ped_cohorts <- master_genomics_hist_file %>%
  pull(cohort) %>%
  unique()

# this is not fully functional
# use pathology diagnosis mapping where short/broad histology is NA 
# source: https://github.com/d3b-center/D3b-codes/tree/7fbb66c8eb1dc9d8086d0c903333810a2bd9d422/OpenPBTA_v20_release_QC/input
# pathology_diagnosis_code <- read.delim(file.path(data_dir, "pathology_diagnosis_for_subtyping.tsv"))
# colnames(pathology_diagnosis_code) <- paste0('map_', colnames(pathology_diagnosis_code))
# master_genomics_hist_file = master_genomics_hist_file %>%
#   left_join(pathology_diagnosis_code, by = c("pathology_diagnosis" = "map_pathology_diagnosis"))
# master_genomics_hist_file <- master_genomics_hist_file %>%
#   mutate(short_histology = ifelse(is.na(short_histology), yes = map_short_histology, no = short_histology),
#          broad_histology = ifelse(is.na(broad_histology), yes = map_broad_histology, no = broad_histology)) 

# # subset to cols of interest
# add_cols <- setdiff(hist_cols_used, colnames(master_genomics_hist_file))
# if(length(add_cols) > 0){
#   master_genomics_hist_file[,add_cols] <- NA
# }
# master_genomics_hist_file <- master_genomics_hist_file %>% 
#   dplyr::select(hist_cols_used)

# create histology subset
ped_hist_file <- create_histologies_subset(hist_file = master_genomics_hist_file, 
                                           group_filter = pediatric_cancer,
                                           cohort_filter = ped_cohorts,
                                           output_dir = pediatric_cancer_dir,
                                           prefix = "pediatric_tumors")

# tpm matrices
tpm_file <- master_genomics_files[grep("tpm", master_genomics_files)]
tpm_mat <- file.path(tpm_file)
tpm_mat <- readRDS(tpm_mat)

# count matrices
counts_file <- master_genomics_files[grep("expected_count", master_genomics_files)]
counts_mat <- file.path(counts_file)
counts_mat <- readRDS(counts_mat)

# pediatric rnaseq matrices
create_rnaseq_matrices(hist_file = ped_hist_file, 
                       tpm_mat = tpm_mat, 
                       counts_mat = counts_mat, 
                       output_dir = pediatric_cancer_dir, 
                       prefix = "pediatric_tumors")

# input files
# copy number
cnv_file <- master_genomics_files[grep("consensus_wgs_plus_cnvkit_wxs", master_genomics_files)]
cnv_data <- data.table::fread(cnv_file)

# mutations
snv_file <- master_genomics_files[grep("snv", master_genomics_files)]
# skip comment lines while reading
snv_data <- data.table::fread(cmd = paste0("gzcat ", snv_file, " | grep -v '^#' ")) 

# pediatric mutations/copy number 
create_non_rnaseq_matrices(hist_file = ped_hist_file,
                           cnv_data = cnv_data,
                           snv_data = snv_data,
                           output_dir = pediatric_cancer_dir,
                           prefix = "pediatric_tumors")

# pediatric fusions
fusion_file <- master_genomics_files[grep("starfusion", master_genomics_files)]
star_fusion <- data.table::fread(fusion_file)
fusion_file <- master_genomics_files[grep("arriba", master_genomics_files)]
arriba_fusion <- data.table::fread(fusion_file)
create_fusion_files(star_fusion = star_fusion, 
                    arriba_fusion = arriba_fusion, 
                    hist_file = ped_hist_file,
                    output_dir = pediatric_cancer_dir,
                    prefix = "pediatric_tumors")


# pediatric tumor vs normal diff expr (for oncogrid)
diff_expr(hist_file = ped_hist_file,
          pediatric_cancer_dir = pediatric_cancer_dir, 
          normal_tissue_dir = normal_tissue_dir, 
          output_dir = pediatric_cancer_dir, 
          prefix = "pediatric_tumors")

# 4. full pediatric data 
# create histology subset
pediatric_cohort_dir <- file.path(data_dir, 'pediatric_data_all')
ped_hist_file <- create_histologies_subset(hist_file = master_genomics_hist_file, 
                                           group_filter = NULL,
                                           cohort_filter = ped_cohorts,
                                           output_dir = pediatric_cohort_dir,
                                           prefix = "pediatric_tumors")

# pediatric rnaseq matrices
create_rnaseq_matrices(hist_file = ped_hist_file, 
                       tpm_mat = tpm_mat, 
                       counts_mat = counts_mat, 
                       output_dir = pediatric_cohort_dir, 
                       prefix = "pediatric_tumors")

# pediatric mutations/copy number 
create_non_rnaseq_matrices(hist_file = ped_hist_file,
                           cnv_data = cnv_data,
                           snv_data = snv_data,
                           output_dir = pediatric_cohort_dir,
                           prefix = "pediatric_tumors")

# pediatric fusions
create_fusion_files(star_fusion = star_fusion, 
                    arriba_fusion = arriba_fusion, 
                    hist_file = ped_hist_file,
                    output_dir = pediatric_cohort_dir,
                    prefix = "pediatric_tumors")
