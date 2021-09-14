# Author: Run Jin
# Classify SNV mutation results into 4 tier [@doi:10.1016/j.jmoldx.2016.10.002]

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(readr)
})

# function to get all directories

# Define Directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
input_dir <- file.path(patient_dir, "output", "oncokb_analysis")

# Read in files necessary for analyses
oncokb_anno <- readr::read_tsv(file.path(input_dir, "oncokb_consensus_annotated.txt"))
cancer_genes <- readRDS(file.path(data_dir, 'cancer_gene_list.rds'))

# find the oncogenes tumor suppressor genes in OMPARE knowledge base
oncogenes_tsg <- cancer_genes %>%
  filter(type %in% c("Is.Oncogene", "Oncogene", "Is.Tumor.Suppressor.Gene", "TumorSuppressorGene")) %>%
  .$Gene_Symbol %>% unique

# calculate VAF for subsequent classification 
oncokb_anno <- oncokb_anno %>% 
  dplyr::mutate(VAF = t_alt_count/(t_alt_count+t_ref_count))

# Filter for Tier I
oncokb_anno_tier1 <- oncokb_anno %>% 
  # filter to variation of interest 
  dplyr::filter(Variant_Classification %in% c('Missense_Mutation', 'Nonsense_Mutation', # filter the files 
                              'Frame_Shift_Del', 'Frame_Shift_Ins', 
                              'In_Frame_Del', 'In_Frame_Ins', 
                              'Splice_Region', 'Splice_Site')) %>%
  # filter based on OncoKB annotation level 
  dplyr::filter(HIGHEST_LEVEL %in% c("LEVEL_1","LEVEL_2","LEVEL_3A", "LEVEL_3B", "LEVEL_R1") | HIGHEST_DX_LEVEL %in% c("LEVEL_Dx1","LEVEL_Dx2") | HIGHEST_PX_LEVEL %in% c("LEVEL_Px1","LEVEL_Px2")) %>%
  # contain only variants that are <1% (less likely to be germline)
  dplyr::filter(gnomAD_AF<0.01| is.na(gnomAD_AF)) %>%
  # filter on VAF over 5% to avoid background noise
  dplyr::filter(VAF>0.05) %>%
  # filter to contain only mutations that has COSMIC ID
  dplyr::mutate(has_COSMIC = dplyr::case_when(
    stringr::str_detect(Existing_variation,
                        "COSM") ~ "Yes",
    TRUE ~ "No"
  )) %>%
  dplyr::filter(has_COSMIC == "Yes") %>%
  # filter to contain only mutations with citations
  dplyr::filter(!is.na(CITATIONS)) %>% 
  # filter to TSG lists in OMPARE database
  dplyr::filter(Hugo_Symbol %in% oncogenes_tsg ) %>% 
  dplyr::select(-has_COSMIC) %>% 
  dplyr::mutate(tier_classification = "TIER I")
  

# Filter for Tier II
oncokb_anno_tier2 <- oncokb_anno %>% 
  # filter to variation of interest 
  dplyr::filter(Variant_Classification %in% c('Missense_Mutation', 'Nonsense_Mutation', # filter the files 
                                              'Frame_Shift_Del', 'Frame_Shift_Ins', 
                                              'In_Frame_Del', 'In_Frame_Ins', 
                                              'Splice_Region', 'Splice_Site')) %>%
  # filter based on OncoKB annotation level 
  dplyr::filter(HIGHEST_LEVEL %in% c("LEVEL_4","LEVEL_R2") | HIGHEST_DX_LEVEL %in% c("LEVEL_Dx3") | HIGHEST_PX_LEVEL %in% c("LEVEL_Px3")) %>%
  # contain only variants that are <1% (less likely to be germline)
  dplyr::filter(gnomAD_AF<0.01| is.na(gnomAD_AF)) %>%
  # filter on VAF over 5% to avoid background noise
  dplyr::filter(VAF>0.05) %>%
  # filter to contain only mutations that has COSMIC ID
  dplyr::mutate(has_COSMIC = dplyr::case_when(
    stringr::str_detect(Existing_variation,
                        "COSM") ~ "Yes",
    TRUE ~ "No"
  )) %>%
  dplyr::filter(has_COSMIC == "Yes") %>%
  # filter to contain only mutations with citations
  dplyr::filter(!is.na(CITATIONS)) %>% 
  # filter to TSG lists in OMPARE database
  dplyr::filter(Hugo_Symbol %in% oncogenes_tsg ) %>% 
  dplyr::select(-has_COSMIC) %>% 
  dplyr::mutate(tier_classification = "TIER II")


# Filter for Tier III
oncokb_anno_tier3 <- oncokb_anno %>% 
  # filter to variation of interest 
  dplyr::filter(Variant_Classification %in% c('Missense_Mutation',  # filter the files 
                                              'In_Frame_Del', 'In_Frame_Ins')) %>%
  # filter to samples that do not have OncoKB annotation
  dplyr::filter(is.na(HIGHEST_LEVEL)) %>%
  dplyr::filter(is.na(HIGHEST_DX_LEVEL)) %>%
  dplyr::filter(is.na(HIGHEST_PX_LEVEL)) %>%
  # contain only variants that are <1% (less likely to be germline)
  dplyr::filter(gnomAD_AF<0.01| is.na(gnomAD_AF)) %>%
  # filter on VAF over 5% to avoid background noise
  dplyr::filter(VAF>0.05) %>%
  dplyr::mutate(tier_classification = "TIER III")

# Filter for Tier IV
oncokb_anno_tier4 <- oncokb_anno %>% 
  # filter to variation of interest 
  dplyr::filter(Variant_Classification =='Missense_Mutation') %>%
  # filter to samples that do not have OncoKB annotation
  dplyr::filter(is.na(HIGHEST_LEVEL)) %>%
  dplyr::filter(is.na(HIGHEST_DX_LEVEL)) %>%
  dplyr::filter(is.na(HIGHEST_PX_LEVEL)) %>%
  # contain only variants that are <1% (less likely to be germline)
  dplyr::filter(gnomAD_AF>=0.01) %>%
  # filter on VAF over 5% to avoid background noise
  dplyr::filter((VAF>=0.4 & VAF<=0.6) | VAF>=0.9) %>%
  dplyr::mutate(tier_classification = "TIER IV")

# combine all the variants in MAF that contain Tier information 
oncokb_anno_tiers <- bind_rows(oncokb_anno_tier1, oncokb_anno_tier2, oncokb_anno_tier3, oncokb_anno_tier4) %>%
  dplyr::select(Hugo_Symbol, Variant_Classification, HGVSp_Short, tier_classification)

# annotate TIER information to Key clinical findings output 

key_clinical_findings_output_modified <- data.frame()

for(i in 1:nrow(oncokb_anno_tiers)) {
  # for each mutation, find the matching entries in the key clinical findings table (if present)
  each <- oncokb_anno_tiers[i,]
  gene_symbol <- each$Hugo_Symbol
  variant_classification <- each$Variant_Classification
  hgvs_short <- each$HGVSp_Short
  tier_class <- each$tier_classification
  
  # make changes to variant property column when TIER information is available
  key_clinical_findings_output_added <- key_clinical_findings_output %>% 
    dplyr::filter(grepl(variant_classification, Details) & grepl(hgvs_short, Details)) %>%
    dplyr::mutate(Variant_Properties = paste0(tier_class, "; ", Variant_Properties))
  
  # save the edited entries
  key_clinical_findings_output_modified <- key_clinical_findings_output_modified %>% rbind(key_clinical_findings_output_added)
}

# merge the modified rows back to the original table 
key_clinical_findings_output_unchanged <- key_clinical_findings_output %>% 
  dplyr::filter(!Details %in% key_clinical_findings_output_modified$Details)

key_clinical_findings_output <- bind_rows(key_clinical_findings_output_modified, key_clinical_findings_output_unchanged)


