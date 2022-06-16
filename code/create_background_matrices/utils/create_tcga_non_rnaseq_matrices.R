# function to create TCGA matrices

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(tidyverse)
  library(reshape2)
  library(dplyr)
})

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
module_dir <- file.path(root_dir, "code", "create_background_matrices")

# source functions
source(file.path(module_dir, "utils", "filter_cnv.R"))
source(file.path(module_dir, "utils", "filter_mutations.R"))

# subset to cancer genes 
cancer_genes <- readRDS(file.path(data_dir, "cancer_gene_list.rds"))

create_tcga_non_rnaseq_matrices <- function(hist_file,
                                            cohort_filter = c("TCGA"), 
                                            group_filter, 
                                            output_dir, 
                                            prefix){
  
  # create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # set project to query
  group_filter <- hist_file %>%
    pull(short_histology) %>%
    unique()
  project <- paste(cohort_filter, group_filter, sep = "-")
  
  # add patient barcode
  hist_file <- hist_file%>%
    dplyr::select(Kids_First_Biospecimen_ID, cohort_participant_id, cohort, sample_id)

  # mutations (hg38)
  query <- GDCquery(
    project = project, 
    data.category = "Simple Nucleotide Variation", 
    access = "open", 
    legacy = FALSE, 
    data.type = "Masked Somatic Mutation", 
    workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
  )
  GDCdownload(query)
  maf <- GDCprepare(query)
  maf <- maf[grep('mutect2', maf$callers),]
  maf$sample_id <- gsub('[D]-[0-9A-Z]{4}-[0-9]{2}', '', maf$Tumor_Sample_Barcode)
  maf <- maf %>%
    mutate(Kids_First_Biospecimen_ID = Tumor_Sample_Barcode) %>%
    inner_join(hist_file %>%
                 dplyr::select(-c(Kids_First_Biospecimen_ID)), by = "sample_id")
  
  # filter mutations
  maf <- filter_mutations(myMutData = maf, myCancerGenes = cancer_genes)
  saveRDS(maf, file = file.path(output_dir, paste(prefix, 'mutation_filtered.rds', sep = "_")))

  # # copy number
  # query <- GDCquery(project = project,
  #                   data.category = "Copy Number Variation",
  #                   data.type = "Gene Level Copy Number Scores",
  #                   legacy = TRUE,
  #                   workflow.type = "ASCAT2")
  # GDCdownload(query)
  # cnv_data <- GDCprepare(query, summarizedExperiment = F)
  # 
  # # convert matrix to long format
  # cnv_data <- cnv_data %>%
  #   dplyr::select(-c(`Gene ID`)) %>% 
  #   dplyr::rename("ensembl" = "Gene Symbol",
  #                 "cytoband" = "Cytoband") %>%
  #   gather('Kids_First_Biospecimen_ID', 'copy_number', -c("ensembl", "cytoband")) 
  # 
  # # format identifiers so we can map to RNA
  # cnv_data$sample_id <- gsub('[D]-[0-9A-Z]{4}-[0-9]{2}', '', cnv_data$Kids_First_Biospecimen_ID)
  # 
  # # only keep gain/loss and filter to samples of interest
  # cnv_data <- cnv_data %>%
  #   filter(copy_number != 0) %>%
  #   inner_join(hist_file %>%
  #                dplyr::select(-c(Kids_First_Biospecimen_ID)), by = "sample_id") %>%
  #   mutate(status = ifelse(copy_number == 1, 'Gain', 'Loss')) %>%
  #   unique()
  # 
  # # map gene symbols
  # cnv_data <- cnv_data %>%
  #   inner_join(annot, by = c('ensembl' = 'original_ensembl_gene_id')) %>%
  #   dplyr::rename("hgnc_symbol" = "external_gene_name") 
  # 
  # # filter to cancer genes (Oncogenes and TSGs only)
  # cnv_data <- filter_cnv(myCNVData = cnv_data, myCancerGenes = cancer_genes)
  # 
  # cnv_data <- cnv_data %>%
  #   mutate(cohort = cohort_filter,
  #          ploidy = NA) %>%
  #   dplyr::select(Kids_First_Biospecimen_ID, status, copy_number, ploidy, ensembl, hgnc_symbol, cytoband, cohort_participant_id, cohort, sample_id) %>%
  #   unique()
  # saveRDS(cnv_data, file = file.path(output_dir, paste(prefix, 'cnv_filtered.rds', sep = "_")))
}
