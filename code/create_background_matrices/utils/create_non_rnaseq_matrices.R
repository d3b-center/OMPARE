# function to create non-rnaseq matrices for mutational analysis
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(dplyr)
})

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
module_dir <- file.path(root_dir, "code", "create_background_matrices")

# source functions
source(file.path(module_dir, "utils", "filter_cnv.R"))
source(file.path(module_dir, "utils", "filter_mutations.R"))

# cancer genes 
cancer_genes <- readRDS(file.path(data_dir, "cancer_gene_list.rds"))

create_non_rnaseq_matrices <- function(hist_file, cnv_data, snv_data, output_dir, prefix){
  
  # create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # subset to id columns
  hist_file <- hist_file %>%
    dplyr::select(Kids_First_Biospecimen_ID, sample_id, cohort, cohort_participant_id) %>%
    unique()
  
  # consensus mutations 
  snv_data <- snv_data %>%
    mutate(Kids_First_Biospecimen_ID = Tumor_Sample_Barcode) %>%
    inner_join(hist_file, by = "Kids_First_Biospecimen_ID")

  # filter mutations
  snv_data <- filter_mutations(myMutData = snv_data,  myCancerGenes = cancer_genes)
  saveRDS(snv_data, file = file.path(output_dir, paste(prefix, 'mutation_filtered.rds', sep = "_")))

  # read cnvkit
  if(!is.null(cnv_data)){
    cnv_data <- cnv_data %>%
      rename("Kids_First_Biospecimen_ID" = "biospecimen_id") %>%
      inner_join(hist_file, by = "Kids_First_Biospecimen_ID") 
    
    # modify copy number status
    cnv_data <- cnv_data %>%
      mutate(status = case_when(status == "gain" ~ "Gain",
                                status == "loss" ~ "Loss",
                                status == "neutral" ~ "Neutral",
                                status == "deep deletion" ~ "Complete Loss",
                                status == "amplification" ~ "Amplification"))
    
    # filter to cancer genes (Oncogenes and TSGs only)
    cnv_data <- filter_cnv(myCNVData = cnv_data %>% dplyr::rename("hgnc_symbol" = "gene_symbol"), 
                           myCancerGenes = cancer_genes)
    cnv_data <- cnv_data %>%
      dplyr::select(Kids_First_Biospecimen_ID, status, copy_number, ploidy, ensembl, hgnc_symbol, cytoband, cohort_participant_id, cohort, sample_id) %>%
      unique()
    saveRDS(cnv_data, file = file.path(output_dir, paste(prefix, 'cnv_filtered.rds', sep = "_")))
  }
}