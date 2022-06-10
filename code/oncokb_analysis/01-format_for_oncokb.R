# filter and format cnv and fusion files for oncokb annotator
suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
})

# arguments
option_list <- list(
  make_option(c("--patient_of_interest"), type = "character",
              help = "cohort participant id for patient of interest"),
  make_option(c("--master_genomics_dir"), type = "character",
              help = "directory with master genomics files"),
  make_option(c("--output_dir"), type = "character",
              help = "output directory")
)
opt <- parse_args(OptionParser(option_list = option_list))
patient_of_interest <- opt$patient_of_interest
master_genomics_dir <- opt$master_genomics_dir
output_dir <- opt$output_dir

# set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# output files
cnv_out <- file.path(output_dir, "oncokb_cnv.txt")
fusion_out <- file.path(output_dir, "oncokb_fusion.txt")

# master genomics file list
master_genomics_files  <- list.files(path = master_genomics_dir, all.files = T, full.names = T)

# sample info for sample of interest
sample_info <- master_genomics_files[grep('hist', master_genomics_files)]
sample_info <- data.table::fread(sample_info)
sample_info <- sample_info %>%
  filter(cohort_participant_id == patient_of_interest)

if(!file.exists(cnv_out)){
  # format cnv for cnaAnnotator.py
  cnv_data <- master_genomics_files[grep('cnv', master_genomics_files)]
  cnv_data <- data.table::fread(cnv_data)
  cnv_data <- cnv_data %>%
    filter(biospecimen_id %in% sample_info$Kids_First_Biospecimen_ID)

  # write out matrix of GISTIC like output
  cnv_data <- cnv_data %>%
    mutate(value = case_when(status == "deep deletion" ~ -2,
                              status == "loss" ~ -1,
                              status == "neutral" ~ 0,
                              status == "gain" ~ 1,
                              status == "amplification" ~ 2)) %>% 
    mutate(`Gene Symbol` = gene_symbol,
           `Locus ID` = 0,
           Cytoband = cytoband)  %>%
    dplyr::select(c(`Gene Symbol`, `Locus ID`, Cytoband, value)) %>%
    dplyr::rename(!!patient_of_interest := value)
  write.table(cnv_data, file = cnv_out, quote = F, sep = "\t", row.names = F)
}

if(!file.exists(fusion_out)){
  # format fusions for FusionAnnotator.py 
  # arriba: filter for high confidence fusions only
  arriba_data <- master_genomics_files[grep('arriba', master_genomics_files)]
  arriba_data <- data.table::fread(arriba_data)
  arriba_data <- arriba_data %>%
    filter(tumor_id %in% sample_info$Kids_First_Biospecimen_ID)
  arriba_data <- arriba_data %>%
    filter(confidence == "high") %>%
    dplyr::mutate(Fusion = paste0(gene1, '-', gene2),
           Tumor_Sample_Barcode = patient_of_interest) %>%
    dplyr::select(Tumor_Sample_Barcode, Fusion)
  
  # star fusion: filter for FFPM>0.1 AND LargeAnchorSupport (YES_LDAS)
  starfusion_data <- master_genomics_files[grep('starfusion', master_genomics_files)]
  starfusion_data <- data.table::fread(starfusion_data)
  starfusion_data <- starfusion_data %>%
    filter(tumor_id %in% sample_info$Kids_First_Biospecimen_ID) %>%
    mutate(Fusion = gsub('--','-',FusionName),
           Tumor_Sample_Barcode = patient_of_interest) %>%
    filter(FFPM > 0.1,
           LargeAnchorSupport == "YES_LDAS") %>%
    dplyr::select(Tumor_Sample_Barcode, Fusion)
  
  # combine both and write out non-redundant arriba + star-fusion calls
  fusion_data <- unique(rbind(starfusion_data, arriba_data))
  write.table(fusion_data, file = fusion_out, quote = F, sep = "\t", row.names = F)
}

  
