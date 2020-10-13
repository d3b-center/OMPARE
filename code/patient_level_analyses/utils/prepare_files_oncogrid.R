# Author: Komal S. Rathi
# Date: 08/10/2020
# Function: Script to generat Oncogrid matrix/additional files using CNV, SNV, Fusion and Expression data
suppressPackageStartupMessages(library(tidyverse))

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# oncogrid directory
oncogrid_path <- file.path(ref_dir, 'oncogrid')
oncogrid_path_input <- file.path(oncogrid_path, 'input')
oncogrid_path_output <- file.path(oncogrid_path, 'output')

# pnoc008 directory
pnoc008_dir <- file.path(ref_dir, 'pnoc008')

## Oncoprint matrix
# cohort 3 matrix
cohort3.matrix <- read.delim(file.path(oncogrid_path_output, 'oncoprint.cohort3.txt'), check.names = F)
colnames(cohort3.matrix)[1] <- 'Sample' 
tmb_mat <- read.delim(file.path(oncogrid_path_input, 'sample_cohort3.TMB'), check.names = F)
tmb_mat <- tmb_mat %>% 
  dplyr::rename(Sample = 1)
annot_info <- read.delim(file.path(oncogrid_path_input, 'sample_cohort3.info'), check.names = F)
annot_info <- annot_info %>% 
  dplyr::rename(Sample = 1)
split_samples <- read.delim(file.path(oncogrid_path_input, 'sample_cohort3_split_samples.txt'), check.names = F)
split_samples <- split_samples %>% 
  dplyr::rename(Sample = 1)

# read reference gene lists
snv <- read.delim(file.path(oncogrid_path_input, 'snv-genes'), header = F)
fusion <- read.delim(file.path(oncogrid_path_input, 'fusion_genes'), header = F)
cnv <- read.delim(file.path(oncogrid_path_input, 'copy_number_gene'), header = F)
deg <- read.delim(file.path(oncogrid_path_input, 'all_cnv_tgen_genes'), header = F)

# fill in details for PNOC008 patients
# 1. get degene info PNOC008 patientss vs GTEx Brain
fname <- file.path(pnoc008_dir, 'pnoc008_vs_gtex_brain_degs.rds')
genes_df <- readRDS(fname)
deg_genes <- genes_df %>%
  mutate(label = ifelse(diff_expr == "up", "OVE", "UNE"),
         Gene_name = genes) %>%
  filter(Gene_name %in% deg$V1) %>%
  dplyr::select(sample_name, Gene_name, label) %>%
  unique()

# 2. get cnv info
fname <- file.path(pnoc008_dir, 'pnoc008_cnv_filtered.rds')
cnv_genes <- readRDS(fname)
cnv_genes <- cnv_genes %>%
  mutate(label = ifelse(Alteration_Type == "Gain", "GAI", "LOS")) %>%
  filter(Gene %in% cnv$V1) %>%
  mutate(Gene_name = Gene,
         sample_name = Kids_First_Biospecimen_ID) %>%
  dplyr::select(sample_name, Gene_name, label) %>%
  unique()

# 3. get snv info
fname <- file.path(pnoc008_dir, 'pnoc008_consensus_mutation_filtered.rds')
mut_genes <- readRDS(fname)
mut_genes <- mut_genes %>%
  filter(!Alteration_Type %in% c("3'Flank", "5'Flank", "3'UTR", "5'UTR", "IGR", "Intron", "RNA")) %>%
  mutate(label = case_when(Alteration_Type %in% "Missense_Mutation" ~ "MIS",
                           Alteration_Type %in% "Nonsense_Mutation" ~ "NOS",
                           Alteration_Type %in% "Frame_Shift_Del" ~ "FSD",
                           Alteration_Type %in% "Frame_Shift_Ins" ~ "FSI",
                           Alteration_Type %in% "In_Frame_Del" ~ "IFD",
                           Alteration_Type %in% "Splice_Site" ~ "SPS")) %>%
  filter(Gene %in% snv$V1) %>%
  mutate(Gene_name = Gene,
         sample_name = Kids_First_Biospecimen_ID) %>%
  dplyr::select(sample_name, Gene_name, label) %>%
  unique()

# 4. get fusion info
fname <- file.path(pnoc008_dir, 'pnoc008_fusions_filtered.rds')
fus_genes <- readRDS(fname)
fus_genes <- fus_genes %>%
  mutate(label = "FUS") %>%
  filter(Gene %in% fusion$V1) %>%
  mutate(Gene_name = Gene,
         sample_name = Kids_First_Biospecimen_ID) %>%
  dplyr::select(sample_name, Gene_name, label) %>%
  unique()

# combine fus + snv
snv_fus <- rbind(mut_genes, fus_genes)

# combine deg + cnv
cnv_deg <- rbind(cnv_genes, deg_genes)

# uniquify rows
snv_fus <- snv_fus %>%
  group_by(sample_name, Gene_name) %>%
  summarise(label = paste0(label, collapse = ';'))
cnv_deg <- cnv_deg %>%
  group_by(sample_name, Gene_name) %>%
  summarise(label = paste0(label, collapse = ';'))

# convert to matrix
snv_fus <- snv_fus %>%
  spread(key = Gene_name, value = 'label') %>%
  column_to_rownames('sample_name')
cnv_deg <- cnv_deg %>%
  spread(key = Gene_name, value = 'label') %>%
  column_to_rownames('sample_name')

# add an * to common genes 
colnames(cnv_deg) <- ifelse(colnames(cnv_deg) %in% colnames(snv_fus), paste0(colnames(cnv_deg),'*'), colnames(cnv_deg))

# merge both matrices
pnoc008_oncogrid_mat <- snv_fus %>%
  rownames_to_column('Sample') %>%
  full_join(cnv_deg %>%
              rownames_to_column('Sample'), by = "Sample")
write.table(pnoc008_oncogrid_mat, file = file.path(oncogrid_path_output, 'oncoprint.pnoc008.txt'), quote = F, sep = "\t", row.names = F)

# add cohort3 matrix to pnoc008 matrix
cohort3.pnoc008.matrix <- plyr::rbind.fill(cohort3.matrix, pnoc008_oncogrid_mat)
write.table(cohort3.pnoc008.matrix, file = file.path(oncogrid_path_output, 'oncoprint.cohort3.pnoc008.txt'), quote = F, sep = "\t", row.names = F)

## TMB
# add TMB info
tmb_mat_p <- readRDS(file.path(pnoc008_dir, 'pnoc008_tmb_scores.rds'))
colnames(tmb_mat_p) <- colnames(tmb_mat)
tmb_mat <- unique(rbind(tmb_mat, tmb_mat_p))
tmb_mat <- tmb_mat[match(cohort3.pnoc008.matrix$Sample, tmb_mat$Sample),]
write.table(tmb_mat, file = file.path(oncogrid_path_output, 'tmb.cohort3.pnoc008.txt'), quote = F, sep = "\t", row.names = F)

## Annotation
# add annotation info
annot_info <- annot_info %>% 
  dplyr::rename(Sample = 1)
annot_info_p <- data.frame(Sample = grep('PNOC', cohort3.pnoc008.matrix$Sample, value = T),
                           Sequencing_Experiment = "WXS,RNA-Seq",
                           Cohort = "PNOC008",
                           Tumor_Descriptor = "Primary",
                           Integrated_Diagnosis = "High_grade_glioma",
                           OS_Status = "LIVING")
annot_info <- unique(rbind(annot_info, annot_info_p))
annot_info <- annot_info[match(cohort3.pnoc008.matrix$Sample, annot_info$Sample),]
write.table(annot_info, file = file.path(oncogrid_path_output, 'annotation.cohort3.pnoc008.txt'), quote = F, sep = "\t", row.names = F)

## Split samples
split_samples_p <- data.frame(Sample = grep('PNOC', cohort3.pnoc008.matrix$Sample, value = T),
                              split_samples = "Others")
split_samples <- unique(rbind(split_samples, split_samples_p))
split_samples <- split_samples[match(cohort3.pnoc008.matrix$Sample, split_samples$Sample),]
write.table(split_samples, file = file.path(oncogrid_path_output, 'split_samples.cohort3.pnoc008.txt'), quote = F, sep = "\t", row.names = F)

