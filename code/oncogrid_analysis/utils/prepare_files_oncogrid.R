# Function: script to generate oncogrid matrix/additional files using CNV, SNV, Fusion and Expression data
suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
})

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")

# oncogrid directory
oncogrid_path <- file.path(data_dir, "oncogrid")

# read reference gene lists (from PNOC003)
snv <- read.delim(file.path(oncogrid_path, "snv-genes.tsv"), header = F)
fusion <- read.delim(file.path(oncogrid_path, "fusion_genes.tsv"), header = F)
cnv <- read.delim(file.path(oncogrid_path, "copy_number_gene.tsv"), header = F)
deg <- read.delim(file.path(oncogrid_path, "all_cnv_tgen_genes.tsv"), header = F)

prepare_files_oncogrid <- function(sample_info, pediatric_cancer_dir, output_dir){
  
  # we will use cohort_participant_id to uniquely identify samples
  # sample of interest
  patient_of_interest <- unique(sample_info$cohort_participant_id)
  
  # get all data from pediatric tumors 
  pediatric_cancer_files <- list.files(path = pediatric_cancer_dir, full.names = T)
  
  # 1. get degene info (pediatric tumors vs normals)
  deg_genes <- pediatric_cancer_files[grep('degs', pediatric_cancer_files)]
  deg_genes <- readRDS(deg_genes)
  deg_genes <- deg_genes %>% 
    ungroup() %>%
    dplyr::mutate(label = ifelse(diff_expr == "up", "OVE", "UNE"),
                  Gene_name = gene_symbol,
                  sample = cohort_participant_id) %>%
    filter(Gene_name %in% deg$V1) %>%
    dplyr::select(sample, Gene_name, label) %>%
    unique()
  
  # 2. get cnv info
  cnv_genes <- pediatric_cancer_files[grep('cnv_filtered', pediatric_cancer_files)]
  cnv_genes <- readRDS(cnv_genes)
  cnv_genes <- cnv_genes %>% 
    ungroup() %>%
    plyr::mutate(label = ifelse(status %in% c("Gain", "Amplification"), "GAI", "LOS")) %>%
    filter(hgnc_symbol %in% cnv$V1) %>%
    dplyr::mutate(Gene_name = hgnc_symbol,
                  sample = cohort_participant_id) %>%
    dplyr::select(sample, Gene_name, label) %>%
    unique()
  
  # 3. get snv info
  mut_genes <- pediatric_cancer_files[grep('mutation_filtered', pediatric_cancer_files)]
  mut_genes <- readRDS(mut_genes)
  mut_genes <- mut_genes %>% 
    ungroup() %>%
    mutate(sample = cohort_participant_id,
           Gene_name = Hugo_Symbol) %>%
    filter(!Variant_Classification %in% c("3'Flank", "5'Flank", "3'UTR", "5'UTR", "IGR", "Intron", "RNA")) %>%
    dplyr::mutate(label = case_when(Variant_Classification %in% "Missense_Mutation" ~ "MIS",
                                    Variant_Classification %in% "Nonsense_Mutation" ~ "NOS",
                                    Variant_Classification %in% "Frame_Shift_Del" ~ "FSD",
                                    Variant_Classification %in% "Frame_Shift_Ins" ~ "FSI",
                                    Variant_Classification %in% "In_Frame_Del" ~ "IFD",
                                    Variant_Classification %in% "Splice_Site" ~ "SPS")) %>%
    filter(Gene_name %in% snv$V1) %>%
    dplyr::select(sample, Gene_name, label) %>%
    unique()
  
  # 4. get fusion info
  fus_genes <- pediatric_cancer_files[grep('fusion_filtered', pediatric_cancer_files)]
  fus_genes <- readRDS(fus_genes)
  fus_genes <- fus_genes %>%
    ungroup() %>%
    separate_rows(fusion_name, sep = "_") %>%
    dplyr::mutate(label = "FUS",
                  Gene_name = fusion_name,
                  sample = cohort_participant_id) %>%
    filter(Gene_name %in% fusion$V1) %>%
    dplyr::select(sample, Gene_name, label) %>%
    unique()
  
  # combine fus + snv
  snv_fus <- rbind(mut_genes, fus_genes)
  
  # combine deg + cnv
  cnv_deg <- rbind(cnv_genes, deg_genes)
  
  # uniquify rows
  snv_fus <- snv_fus %>%
    group_by(sample, Gene_name) %>%
    dplyr::summarise(label = paste0(label, collapse = ';'))
  cnv_deg <- cnv_deg %>%
    group_by(sample, Gene_name) %>%
    dplyr::summarise(label = paste0(label, collapse = ';'))
  
  # convert to matrix
  snv_fus <- snv_fus %>%
    spread(key = Gene_name, value = 'label')%>%
    column_to_rownames('sample')
  cnv_deg <- cnv_deg %>%
    spread(key = Gene_name, value = 'label') %>%
    column_to_rownames('sample')
  
  # add an * to common genes 
  colnames(cnv_deg) <- ifelse(colnames(cnv_deg) %in% colnames(snv_fus), paste0(colnames(cnv_deg),'*'), colnames(cnv_deg))
  
  # merge both matrices
  oncogrid_mat <- snv_fus %>%
    rownames_to_column('Sample') %>%
    full_join(cnv_deg %>%
                rownames_to_column('Sample'), by = "Sample")
  saveRDS(oncogrid_mat, file = file.path(output_dir, "oncogrid_input_matrix.rds"))
  
  # TMB info - TBD
  
  # add annotation info
  annot_info <- pediatric_cancer_files[grep('histologies', pediatric_cancer_files)]
  annot_info <- read.delim(annot_info)
  annot_info <- annot_info %>%
    group_by(cohort_participant_id) %>%
    dplyr::summarise(Sequencing_Experiment = toString(sort(unique(experimental_strategy))),
                     Sample = cohort_participant_id,
                     Cohort = cohort,
                     Tumor_Descriptor = tumor_descriptor,
                     Short_Histology = short_histology,
                     OS_Status = OS_status) %>%
    unique()
  annot_info <- annot_info %>%
    filter(Sample %in% oncogrid_mat$Sample) %>%
    mutate(Cohort = ifelse(cohort_participant_id == patient_of_interest, 'POI', Cohort))
  annot_info <- annot_info[match(oncogrid_mat$Sample, annot_info$Sample),]
  annot_info$Tumor_Descriptor[annot_info$Tumor_Descriptor == "Primary Tumor"] <- "Primary"
  annot_info$OS_Status[annot_info$OS_Status == "Unknown"] <- NA
  saveRDS(annot_info, file = file.path(output_dir, "oncogrid_input_annotation.rds"))
}
