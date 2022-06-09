# Author: Komal S. Rathi
# Function: for all PNOC008 data: 
# create a matrix of expression
# create metadata using clinical files
# create full summary data files of cnv, mutations and fusions
# this is to be run everytime a new patient comes in - before generating the report
# these summary files are required by mutational analysis and oncogrid

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(stringr)
  library(tidyverse)
  library(GenomicRanges)
  library(optparse)
})

option_list <- list(
  make_option(c("--clin_file"), type = "character",
              default = NULL,
              help = "Manifest file (.xlsx)")
)
# parameters to pass
opt <- parse_args(OptionParser(option_list = option_list))
clin_file <- opt$clin_file

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
results_dir <- file.path(root_dir, "results") # contains all patient output directories

# source functions
source(file.path(root_dir, "code", "utils", "filter_fusions.R"))
source(file.path(root_dir, "code", "utils", 'filter_mutations.R'))
source(file.path(root_dir, "code", "utils", 'filter_cnv.R'))
source(file.path(root_dir, "code", "utils", 'collapse_rnaseq.R'))

# create output directory
pnoc008_dir <- file.path(data_dir, "pnoc008")

# reference data
# cancer genes 
cancer_genes <- readRDS(file.path(data_dir, 'cancer_gene_list.rds'))
gene_list <- unique(cancer_genes$Gene_Symbol)

# gencode reference
gencode_v27 <- read.delim(file.path(data_dir, "pnoc008", "gencode.v27.primary_assembly.annotation.txt"))
gencode_v27_pc <- gencode_v27 %>%
  filter(biotype == "protein_coding")

# 1. clinical file (manifest)
pnoc008_clinical <- readxl::read_xlsx(clin_file, sheet = 1)
colnames(pnoc008_clinical) <- gsub('[() ]', '.', colnames(pnoc008_clinical))
pnoc008_clinical <- pnoc008_clinical %>%
  filter_all(any_vars(!is.na(.))) %>%
  mutate(PNOC.Subject.ID = gsub('P-','PNOC008-', PNOC.Subject.ID))
pnoc008_clinical <- pnoc008_clinical %>%
  mutate(subject_id = PNOC.Subject.ID,
         ethnicity = Ethnicity,
         OS_days = Age.at.Diagnosis..in.days.,
         OS_status = toupper(Last.Known.Status),
         sex = Gender,
         cohort_participant_id = Research.ID) %>%
  dplyr::select(subject_id, cohort_participant_id, short_histology, broad_histology, ethnicity, sex, OS_days, OS_status) %>%
  as.data.frame()

# add other identifiers from PBTA base histology
master_genomics <- readr::read_tsv(list.files(path = file.path(data_dir, "master_genomics"), pattern = 'hist', full.names = T)) # data assembly histology file up to current poi
master_genomics <- master_genomics %>%
  mutate(RNA_library = ifelse(Kids_First_Biospecimen_ID == "BS_83KGBX2K", "stranded", RNA_library)) %>% # this sample does not have RNA_library info
  filter(cohort_participant_id %in% pnoc008_clinical$cohort_participant_id) %>%
  dplyr::select(Kids_First_Biospecimen_ID, experimental_strategy, RNA_library, sample_id, cohort, cohort_participant_id, pathology_diagnosis, integrated_diagnosis, molecular_subtype)

# combine both
pnoc008_clinical <- pnoc008_clinical %>%
  left_join(master_genomics, by = 'cohort_participant_id') %>%
  filter(!is.na(Kids_First_Biospecimen_ID),
         Kids_First_Biospecimen_ID != "BS_862NMAR7") %>% # remove duplicated PNOC008-05
  mutate(tmp = Kids_First_Biospecimen_ID) %>%
  column_to_rownames('tmp')

# 2. TPM
pnoc008_tpm <- readRDS(list.files(path = file.path(data_dir, "master_genomics"), pattern = "tpm", recursive = TRUE, full.names = T))
samples_to_use <- intersect(colnames(pnoc008_tpm), pnoc008_clinical$Kids_First_Biospecimen_ID)
pnoc008_tpm <- pnoc008_tpm %>%
  dplyr::select(gene_id, samples_to_use)
pnoc008_tpm <- collapse_rnaseq(pnoc008_tpm)

# filter to protein coding genes
pnoc008_tpm <- pnoc008_tpm %>%
  rownames_to_column('gene_symbol') %>%
  filter(gene_symbol %in% gencode_v27_pc$gene_symbol,
         !grepl('^HIST', gene_symbol)) %>%
  column_to_rownames('gene_symbol')

# counts matrix
pnoc008_counts <- readRDS(list.files(path = file.path(data_dir, "master_genomics"), pattern = "count", recursive = TRUE, full.names = T))
pnoc008_counts <- pnoc008_counts %>%
  dplyr::select(gene_id, samples_to_use)
pnoc008_counts <- collapse_rnaseq(pnoc008_counts)

# filter to protein coding genes
pnoc008_counts <- pnoc008_counts %>%
  rownames_to_column('gene_symbol') %>%
  filter(gene_symbol %in% gencode_v27_pc$gene_symbol,
         !grepl('^HIST', gene_symbol)) %>%
  column_to_rownames('gene_symbol')

# save expression and clinical
saveRDS(pnoc008_tpm, file = file.path(pnoc008_dir, 'pnoc008_tpm_matrix.rds'))
saveRDS(pnoc008_counts, file = file.path(pnoc008_dir, 'pnoc008_counts_matrix.rds'))
readr::write_tsv(pnoc008_clinical, file = file.path(pnoc008_dir, "pnoc008_clinical.tsv"))

# now subset to id columns
pnoc008_clinical_subset <- pnoc008_clinical %>%
  dplyr::select(subject_id, Kids_First_Biospecimen_ID, cohort_participant_id, sample_id, cohort)

# copy number
cnv_data <- list.files(path = file.path(data_dir, "master_genomics"), pattern = "cnv", recursive = TRUE, full.names = T)
cnv_data <- data.table::fread(cnv_data)
cnv_data <- cnv_data %>%
  dplyr::rename("Kids_First_Biospecimen_ID" = "biospecimen_id") %>%
  inner_join(pnoc008_clinical_subset, by = "Kids_First_Biospecimen_ID") 

# modify copy number status
cnv_data <- cnv_data %>%
  mutate(status = case_when(status == "gain" ~ "Gain",
                            status == "loss" ~ "Loss",
                            status == "neutral" ~ "Neutral",
                            status == "deep deletion" ~ "Complete Loss",
                            status == "amplification" ~ "Amplification"))

# filter to cancer genes (Oncogenes and TSGs only)
cnv_data <- filter_cnv(myCNVData = cnv_data %>% 
                         dplyr::rename("hgnc_symbol" = "gene_symbol"), 
                       myCancerGenes = cancer_genes)
cnv_data <- cnv_data %>%
  dplyr::select(subject_id, Kids_First_Biospecimen_ID, status, copy_number, ploidy, ensembl, hgnc_symbol, cytoband, cohort_participant_id, cohort, sample_id) %>%
  unique()
saveRDS(cnv_data, file = file.path(pnoc008_dir, "pnoc008_cnv_filtered.rds"))

# mutations
pnoc008_mutations <- data.table::fread(list.files(path = file.path(data_dir, "master_genomics"), pattern = "snv", recursive = TRUE, full.names = T), skip = 1)
pnoc008_mutations <-  pnoc008_mutations %>%
  mutate(Kids_First_Biospecimen_ID = Tumor_Sample_Barcode) %>%
  inner_join(pnoc008_clinical_subset, by = "Kids_First_Biospecimen_ID")
saveRDS(pnoc008_mutations, file = file.path(pnoc008_dir, "pnoc008_consensus_mutation.rds"))

# filter mutations
pnoc008_mutations <- filter_mutations(myMutData = pnoc008_mutations, myCancerGenes = cancer_genes)
saveRDS(pnoc008_mutations, file = file.path(pnoc008_dir, "pnoc008_consensus_mutation_filtered.rds"))

# fusions
arriba_fusion <- data.table::fread(list.files(path = file.path(data_dir, "master_genomics"), pattern = "arriba", recursive = TRUE, full.names = T))
star_fusion <- data.table::fread(list.files(path = file.path(data_dir, "master_genomics"), pattern = "star", recursive = TRUE, full.names = T))

# filter fusions
star_fusion <- filter_fusions(fusion_data = star_fusion, myCancerGenes = cancer_genes, method = "star-fusion")
star_fusion$method <- "star-fusion"
arriba_fusion <- filter_fusions(fusion_data = arriba_fusion, myCancerGenes = cancer_genes, method = "arriba")
arriba_fusion$method <- "arriba_fusion"
fusion_data <- rbind(star_fusion, arriba_fusion)

# join with histology
pnoc008_fusions <- fusion_data %>%
  dplyr::rename("Kids_First_Biospecimen_ID" = "tumor_id") %>%
  inner_join(pnoc008_clinical_subset, by = "Kids_First_Biospecimen_ID") %>%
  group_by(fusion_name, gene1, gene2, annots, subject_id, Kids_First_Biospecimen_ID, sample_id, cohort, cohort_participant_id) %>%
  summarise(method = toString(method),
            type = toString(type)) %>%
  ungroup()
saveRDS(pnoc008_fusions, file = file.path(pnoc008_dir, "pnoc008_fusions_filtered.rds"))

# cohort level tmb scores
tmb_bed_file <- data.table::fread(file.path(data_dir, 'tmb', 'ashion_confidential_exome_v2_2nt_pad.Gh38.bed'))
colnames(tmb_bed_file)  <- c('chr', 'start', 'end')
tmb_bed_file$chr <- paste0("chr", tmb_bed_file$chr)

# read mutect2 for TMB profile
# for this we have to get all mutect2 files for patients under data/mutect2 and merge them
# because there is no data assembly file for mutect2
pnoc008_mutect2 <- list.files(path = file.path(data_dir, "mutect2"), pattern = 'maf', recursive = TRUE, full.names = T)

# function to merge mutations
colnames_protected <- grep('vep', pnoc008_mutect2, value = T, invert = T) %>% 
  head(1) %>% 
  data.table::fread() %>%
  colnames()
colnames_vep <- grep('vep', pnoc008_mutect2, value = T) %>% 
  head(1) %>% 
  data.table::fread() %>%
  colnames()
common_columns <- intersect(colnames_protected, colnames_vep)
merge_mutations <- function(nm, common_columns){
  x <- data.table::fread(nm)
  x <- x %>% 
    dplyr::select(common_columns)
  if(nrow(x) > 1){
    x <- as.data.frame(x)
    return(x)
  } 
}
pnoc008_mutect2 <- lapply(pnoc008_mutect2, FUN = function(x) merge_mutations(x, common_columns = common_columns))
pnoc008_mutect2 <- data.table::rbindlist(pnoc008_mutect2, fill = T)
pnoc008_mutect2 <- pnoc008_mutect2 %>%
  filter(Tumor_Sample_Barcode %in% pnoc008_clinical_subset$Kids_First_Biospecimen_ID)

# calculate the length of the bed file
bed_length <- 0
for (i in 1:nrow(tmb_bed_file)) {
  distance <- as.numeric(tmb_bed_file[i,3]) - as.numeric(tmb_bed_file[i,2])
  bed_length <- bed_length + distance
}

# mutect2 nonsense and missense mutations
var_class = c('Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins',  'In_Frame_Del', 'In_Frame_Ins')
vaf_cutoff = 0.05
var_count = 3
tumor_depth = 25
pnoc008_mutect2 <- pnoc008_mutect2 %>%
  mutate(vaf = t_alt_count/(t_alt_count+t_ref_count)) %>%
  filter(Variant_Classification %in% var_class,
         t_depth >= tumor_depth,
         vaf >= vaf_cutoff,
         t_alt_count >= var_count) %>%
  dplyr::select(Hugo_Symbol, Variant_Classification, Chromosome, Start_Position, End_Position, Tumor_Sample_Barcode) %>%
  unique()

# intersect with bed file
subject <- with(tmb_bed_file, GRanges(chr, IRanges(start = start, end = end)))
query <- with(pnoc008_mutect2, GRanges(Chromosome, IRanges(start = Start_Position, end = End_Position, names = Hugo_Symbol)))
pnoc008_tmb <- findOverlaps(query = query, subject = subject, type = "within")
pnoc008_tmb <- data.frame(pnoc008_mutect2[queryHits(pnoc008_tmb),], tmb_bed_file[subjectHits(pnoc008_tmb),])

# resolve multiple entries
pnoc008_tmb <- pnoc008_tmb %>%
  dplyr::rename("Kids_First_Biospecimen_ID" = "Tumor_Sample_Barcode") %>%
  inner_join(pnoc008_clinical_subset, by = c("Kids_First_Biospecimen_ID"))

# return the number of filtered variants overlapping with the bed file
pnoc008_tmb <- pnoc008_tmb %>%
  group_by(subject_id) %>%
  mutate(num.mis.non = n()) %>%
  dplyr::select(subject_id, cohort_participant_id, cohort, sample_id, num.mis.non)  %>%
  unique() %>%
  mutate(tmb = num.mis.non*1000000/bed_length) %>%
  unique() %>%
  dplyr::select(-c(num.mis.non))
saveRDS(pnoc008_tmb, file = file.path(pnoc008_dir, "pnoc008_tmb_scores.rds"))

