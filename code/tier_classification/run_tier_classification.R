# Author: Run Jin
# Classify SNV mutation results into 4 tier [@doi:10.1016/j.jmoldx.2016.10.002]

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# Define Directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
module_dir <- file.path(root_dir, "code", "tier_classification")

# source function
source(file.path(module_dir, "utils", "tier_classification.R"))

# Read in files necessary for analyses
input_dir <- file.path(patient_dir, "output", "oncokb_analysis")
if(snv_caller != "all"){
  oncokb_anno <- readr::read_tsv(file.path(input_dir, paste0("oncokb_", snv_caller, "_annotated.txt")))
} else {
  oncokb_anno <- readr::read_tsv(file.path(input_dir, paste0("oncokb_consensus_annotated.txt")))
}

# Read in cancer genes list from OMPARE knowledge base
cancer_genes <- readRDS(file.path(data_dir, 'cancer_gene_list.rds'))

# Read in hotspot database
hotspot_indel <- readr::read_tsv(file.path(data_dir, 'hotspots_database', 'hotspot_database_2017_indel.tsv')) %>% 
  distinct()
hotspot_snv <- readr::read_tsv(file.path(data_dir, 'hotspots_database', 'hotspot_database_2017_snv.tsv')) %>%
  distinct()

# read in COSMIC resistance marker df
cosmic_resistance <- readr::read_tsv(file.path(data_dir, 'CosmicResistanceMutations.tsv')) %>% 
  dplyr::select(`Gene Name`, `AA Mutation`, `CDS Mutation`, `Tier`) %>% 
  dplyr::mutate(`Gene Name` = gsub("\\_.*", "", `Gene Name`)) %>% 
  dplyr::mutate(cosmic_tier_anno = case_when(
    Tier == "1" ~ "Yes",
    TRUE ~ "No"
  )) %>% 
  dplyr::rename(HGVSp_Short = `AA Mutation`, 
                HGVSc = `CDS Mutation`, 
                Hugo_Symbol = `Gene Name`) %>%
  distinct()

# annotate tier classification
if(nrow(oncokb_anno) > 0){
  all_findings_output <- tier_classification(all_findings_output = all_findings_output, 
                                             oncokb_anno = oncokb_anno, 
                                             cancer_genes = cancer_genes, 
                                             hotspot_indel = hotspot_indel, 
                                             hotspot_snv = hotspot_snv,
                                             cosmic_resistance = cosmic_resistance)
}
