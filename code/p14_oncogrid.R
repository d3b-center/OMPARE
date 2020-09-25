# Author: Komal S. Rathi
# Date: 08/10/2020
# Function: Script to generat Oncogrid matrix/additional files using CNV, SNV, Fusion and Expression data
library(tidyverse)

# Directories
# oncogrid directories
oncogrid.path <- file.path('data', 'Reference', 'oncogrid')
oncogrid.path.input <- file.path(oncogrid.path, 'input')
oncogrid.path.output <- file.path(oncogrid.path, 'output')
# pnoc008 directory
pnoc008.path <- file.path('data', 'Reference', 'PNOC008')

## Oncoprint matrix
# cohort 3 matrix
cohort3.matrix <- read.delim(file.path(oncogrid.path.output, 'oncoprint.cohort3.txt'), check.names = F)
colnames(cohort3.matrix)[1] <- 'Sample' 
tmb_mat <- read.delim(file.path(oncogrid.path.input, 'sample_cohort3.TMB'), check.names = F)
tmb_mat <- tmb_mat %>% 
  dplyr::rename(Sample = 1)
annot_info <- read.delim(file.path(oncogrid.path.input, 'sample_cohort3.info'), check.names = F)
annot_info <- annot_info %>% 
  dplyr::rename(Sample = 1)
split_samples <- read.delim(file.path(oncogrid.path.input, 'sample_cohort3_split_samples.txt'), check.names = F)
split_samples <- split_samples %>% 
  dplyr::rename(Sample = 1)

# read reference gene lists
snv <- read.delim(file.path(oncogrid.path.input, 'snv-genes'), header = F)
fusion <- read.delim(file.path(oncogrid.path.input, 'fusion_genes'), header = F)
cnv <- read.delim(file.path(oncogrid.path.input, 'copy_number_gene'), header = F)
deg <- read.delim(file.path(oncogrid.path.input, 'all_cnv_tgen_genes'), header = F)

# fill in details for PNOC008 patients
# 1. get degene info PNOC008 patientss vs GTEx Brain
fname <- file.path(pnoc008.path, 'PNOC008_deg_GTExBrain.rds')
genes.df <- readRDS(fname)
deg.genes <- genes.df %>%
  mutate(label = ifelse(DE == "Up", "OVE", "UNE")) %>%
  filter(Gene_name %in% deg$V1) %>%
  dplyr::select(sample_name, Gene_name, label) %>%
  unique()

# 2. get cnv info from cnvGenes
fname <- file.path(pnoc008.path, 'PNOC008_cnvData_filtered.rds')
cnv.genes <- readRDS(fname)
cnv.genes <- cnv.genes %>%
  mutate(label = ifelse(Alteration_Type == "Gain", "GAI", "LOS")) %>%
  filter(Gene %in% cnv$V1) %>%
  mutate(Gene_name = Gene,
         sample_name = Kids_First_Biospecimen_ID) %>%
  dplyr::select(sample_name, Gene_name, label) %>%
  unique()

# 3. get snv info
fname <- file.path(pnoc008.path, 'PNOC008_consensus_mutData_filtered.rds')
mut.genes <- readRDS(fname)
mut.genes <- mut.genes %>%
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
fname <- file.path(pnoc008.path, 'PNOC008_fusData_filtered.rds')
fus.genes <- readRDS(fname)
fus.genes <- fus.genes %>%
  mutate(label = "FUS") %>%
  filter(Gene %in% fusion$V1) %>%
  mutate(Gene_name = Gene,
         sample_name = Kids_First_Biospecimen_ID) %>%
  dplyr::select(sample_name, Gene_name, label) %>%
  unique()

# combine fus + snv
snv_fus <- rbind(mut.genes, fus.genes)

# combine deg + cnv
cnv_deg <- rbind(cnv.genes, deg.genes)

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
pnoc008.oncogrid.mat <- snv_fus %>%
  rownames_to_column('Sample') %>%
  full_join(cnv_deg %>%
              rownames_to_column('Sample'), by = "Sample")
write.table(pnoc008.oncogrid.mat, file = file.path(oncogrid.path.output, 'oncoprint.pnoc008.txt'), quote = F, sep = "\t", row.names = F)

# add cohort3 matrix to pnoc008 matrix
cohort3.pnoc008.matrix <- plyr::rbind.fill(cohort3.matrix, pnoc008.oncogrid.mat)
write.table(cohort3.pnoc008.matrix, file = file.path(oncogrid.path.output, 'oncoprint.cohort3.pnoc008.txt'), quote = F, sep = "\t", row.names = F)

## TMB
# add TMB info
tmb_mat_p <- readRDS(file.path(pnoc008.path, 'PNOC008_TMBscores.rds'))
colnames(tmb_mat_p) <- colnames(tmb_mat)
tmb_mat <- unique(rbind(tmb_mat, tmb_mat_p))
tmb_mat <- tmb_mat[match(cohort3.pnoc008.matrix$Sample, tmb_mat$Sample),]
write.table(tmb_mat, file = file.path(oncogrid.path.output, 'tmb.cohort3.pnoc008.txt'), quote = F, sep = "\t", row.names = F)

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
write.table(annot_info, file = file.path(oncogrid.path.output, 'annotation.cohort3.pnoc008.txt'), quote = F, sep = "\t", row.names = F)

## Split samples
split_samples_p <- data.frame(Sample = grep('PNOC', cohort3.pnoc008.matrix$Sample, value = T),
                              split_samples = "Others")
split_samples <- unique(rbind(split_samples, split_samples_p))
split_samples <- split_samples[match(cohort3.pnoc008.matrix$Sample, split_samples$Sample),]
write.table(split_samples, file = file.path(oncogrid.path.output, 'split_samples.cohort3.pnoc008.txt'), quote = F, sep = "\t", row.names = F)

