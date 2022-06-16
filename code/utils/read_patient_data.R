# Author: Komal Rathi
# read patient-specific data and save all data to global env

# libraries
suppressPackageStartupMessages({
  library(dplyr)
})

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")

# cancer genes for filtering mutation data
cancer_genes <- readRDS(file.path(data_dir, "cancer_gene_list.rds"))
source('code/utils/filter_mutations.R')

read_patient_data <- function(pediatric_cancer_dir = file.path(data_dir, "pnoc008"), patient_of_interest, mut_only = FALSE, rnaseq_only = FALSE){
  
  # patient dir
  patient_dir <- file.path('results', patient_of_interest)
  
  # patient sample info
  sample_info <- list.files(path = pediatric_cancer_dir, pattern = "tsv", full.names = T)
  sample_info <- data.table::fread(sample_info)
  sample_info <- sample_info %>%
    filter(cohort_participant_id == patient_of_interest)
  readr::write_tsv(x = sample_info, file = file.path(root_dir, 'results', patient_of_interest, 'output', 'sample_info.tsv'))
  assign("sample_info", sample_info, envir = globalenv())
  
  # filtered mutations (use only consensus)
  if(!rnaseq_only){
    maf_file <- list.files(path = file.path(patient_dir, 'simple-variants'), pattern = "consensus", recursive = TRUE, full.names = T)
    full_maf <- data.table::fread(input = maf_file, skip = 1)
    filtered_maf <- filter_mutations(myMutData = full_maf, myCancerGenes = cancer_genes)
  } else {
    filtered_maf <- NULL
  }
  assign("filtered_maf", filtered_maf, envir = globalenv())
  
  # filtered copy number data
  if(!mut_only & !rnaseq_only){
    cnv_file <- list.files(path = pediatric_cancer_dir, pattern = "cnv", full.names = T)
    filtered_cnv <- readRDS(cnv_file)
    filtered_cnv <- filtered_cnv %>%
      filter(Kids_First_Biospecimen_ID %in% sample_info$Kids_First_Biospecimen_ID)
  } else {
    filtered_cnv <- NULL
  }
  assign("filtered_cnv", filtered_cnv, envir = globalenv())
  
  if(!mut_only){
    patient_rna_kfid <- sample_info %>%
      filter(experimental_strategy == "RNA-Seq") %>%
      .$Kids_First_Biospecimen_ID
    
    # filtered fusion data (from both star and arriba)
    fusion_file <- list.files(path = pediatric_cancer_dir, pattern = "fusion_filtered.rds", full.names = T)
    filtered_fusions <- readRDS(fusion_file)
    filtered_fusions <- filtered_fusions %>%
      filter(Kids_First_Biospecimen_ID %in% patient_rna_kfid)
    assign("filtered_fusions", filtered_fusions, envir = globalenv())
    
    # expression data: TPM 
    tpm_file <- list.files(path = pediatric_cancer_dir, pattern = "tpm", full.names = T)
    tpm_data <- readRDS(tpm_file)
    tpm_data <- tpm_data %>%
      select(patient_rna_kfid)
    assign("tpm_data", tpm_data, envir = globalenv())
    
    # expression data: counts 
    counts_file <- list.files(path = pediatric_cancer_dir, pattern = "count", full.names = T)
    count_data <- readRDS(counts_file)
    count_data <- count_data %>%
      select(patient_rna_kfid) 
    assign("count_data", count_data, envir = globalenv()) 
  } else {
    filtered_fusions <- NULL
    tpm_data <- NULL
    count_data <- NULL
  }
  assign("filtered_fusions", filtered_fusions, envir = globalenv())
  assign("tpm_data", tpm_data, envir = globalenv())
  assign("count_data", count_data, envir = globalenv()) 
}

