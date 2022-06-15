# Author: Run Jin
# Generate binary matrix for all samples in the cohort indicating presence/absence
# of mutations in gene list of interest 
# Adapted from https://github.com/d3b-center/d3b-pnoc003-HGG-DMG-omics/
# blob/master/analyses/SNF_mutcooc_and_C2/code/make-matrix.R

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(readr)
  library(reshape2)
})

#### Parse command line options ------------------------------------------------
option_list <- list(
  make_option(c("--ref_cancer_dir"),type="character",
              help="directory with cancer type specific reference tumor subset"),
  make_option(c("--genes"), type = "character",
              help = "file containing the genes of interest for building the matrix (.tsv)"),
  make_option(c("--mat_file"), type = "character",
              help = "binary matrix file with gene vs. sample id indicating mutational status (.tsv)"),
  make_option(c("--sample_id_interest"), type = "character",
              help = "`sample_id` of patient of interest"),
  make_option(c("--analysis_type"), type = "character",
              help = "adult or pediatric")
)
opt <- parse_args(OptionParser(option_list = option_list))
ref_cancer_dir <- opt$ref_cancer_dir
genes <- opt$genes
mat_file <- opt$mat_file
sample_id_interest <- opt$sample_id_interest
analysis_type <- opt$analysis_type

# Read in files necessary for analyses
# for now, we will assume that POI is present in the pediatric subset mutation data only (DGD non-RNA seq cases)
# so this script should work fine as is

# histology file
histology_df <- list.files(path = ref_cancer_dir, pattern = "histologies", full.names = T)
histology_df <- readr::read_tsv(histology_df, guess_max=100000)

# gene file containing all the gene names
gene_df <- readr::read_tsv(file.path(genes))
genes <- unique(gene_df$gene)
onco_genes <- gene_df %>% 
  filter(Is_Oncogene == 'Yes') %>%
  pull(gene)
tsg_genes <- gene_df %>% 
  filter(Is_Tumor_Suppressor_Gene == 'Yes') %>%
  pull(gene)

# Note: currently all subset files are being written as rds
# these files have already been modified to add sample_id and other identifiers using the create_background_matrices module
# in the near future we can revert to the original file extensions

# MAF file
maf_file <- list.files(path = ref_cancer_dir, pattern = "mutation", full.names = T)
maf_file <- readRDS(maf_file)

# Fusion file (currently we don't have fusions for adult data so keep it optional)
fusion_file <- list.files(path = ref_cancer_dir, pattern = "fusion", full.names = T)
if(length(fusion_file) > 0){
  fusion_file <- readRDS(fusion_file) 
}

# SV file (currently we don't have a merged file so let's keep this optional)
# also not available for adult tumors
sv_file <- list.files(path = ref_cancer_dir, pattern = "sv-manta", full.names = T)
if(length(sv_file) > 0){
  sv_file <- readRDS(sv_file) 
}

# CNV file 
cnv_file <- list.files(path = ref_cancer_dir, pattern = "cnv", full.names = T)
cnv_file <- readRDS(cnv_file) 

if(analysis_type == "pediatric"){
  # generate a dataframe with matching RNA and DNA samples 
  # (matched on sample id and aliquot id)
  rna_aliquot_df <- histology_df %>% 
    dplyr::filter(experimental_strategy == "RNA-Seq") %>% 
    dplyr::select(Kids_First_Biospecimen_ID, sample_id, aliquot_id) %>% 
    dplyr::mutate(match_aliquot_id = gsub("\\..*", "", aliquot_id)) %>% 
    dplyr::select(-aliquot_id) %>%
    dplyr::rename(Kids_First_Biospecimen_ID.rna=Kids_First_Biospecimen_ID)
  
  # only keep entries where DNA-RNA match have the same parent aliquot ID or is NA
  histology_df_matched <- histology_df %>% 
    dplyr::filter(experimental_strategy != "RNA-Seq",
                  !is.na(aliquot_id)) %>%
    dplyr::mutate(match_aliquot_id = gsub("\\..*", "", aliquot_id)) %>% 
    dplyr::left_join(rna_aliquot_df) %>%
    dplyr::select(-match_aliquot_id) %>%
    # dplyr::select(Kids_First_Biospecimen_ID, Kids_First_Biospecimen_ID.rna, colnames(histology_df)[2:46]) %>%
    dplyr::rename(Kids_First_Biospecimen_ID.dna =Kids_First_Biospecimen_ID)
  
  # additionally, for each sample ID we need to take independent sample for each sample_id
  # implement independent sample list module in OpenPedCan
  set.seed(2021)
  # first sample them
  histology_df_matched <- histology_df_matched[sample(nrow(histology_df_matched)), ]
  # then take distinct
  histology_df_matched <- histology_df_matched %>%
    dplyr::distinct(sample_id, .keep_all = TRUE) 
  histology_df_matched <- histology_df_matched %>%
    dplyr::select(Kids_First_Biospecimen_ID.dna, Kids_First_Biospecimen_ID.rna, sample_id)
  
  # get samples
  bs_ids <- c(histology_df_matched$Kids_First_Biospecimen_ID.dna,
              histology_df_matched$Kids_First_Biospecimen_ID.rna)
  samples <- histology_df_matched$sample_id
  
  # subset histology to match BS ids
  histology_df <- histology_df %>%
    filter(Kids_First_Biospecimen_ID %in% bs_ids,
           sample_id %in% samples) %>%
    dplyr::select(Kids_First_Biospecimen_ID, sample_id) %>%
    unique()
} else if(analysis_type == "adult"){
  # get ids from cnv
  dna_ids <- cnv_file %>%
    dplyr::select(Kids_First_Biospecimen_ID, sample_id) %>%
    dplyr::rename('Kids_First_Biospecimen_ID.dna' = 'Kids_First_Biospecimen_ID') %>%
    unique()
  
  # get ids from maf
  dna_ids <- dna_ids %>%
    rbind(maf_file %>% 
            dplyr::select(Kids_First_Biospecimen_ID, sample_id) %>%
            dplyr::rename('Kids_First_Biospecimen_ID.dna' = 'Kids_First_Biospecimen_ID')) %>%
    unique() 
  
  # get ids from histology file
  rna_ids <- histology_df %>%
    dplyr::select(Kids_First_Biospecimen_ID, sample_id) %>%
    dplyr::rename('Kids_First_Biospecimen_ID.rna' = 'Kids_First_Biospecimen_ID') %>%
    unique()
  
  # combine all
  histology_df_matched <- dna_ids %>%
    dplyr::left_join(rna_ids) %>%
    unique()
  
  # get samples
  bs_ids <- c(histology_df_matched$Kids_First_Biospecimen_ID.dna,
              histology_df_matched$Kids_First_Biospecimen_ID.rna)
  samples <- histology_df_matched$sample_id
  
  # subset histology to match BS ids
  histology_df <- dna_ids %>%
    dplyr::rename('Kids_First_Biospecimen_ID' = 'Kids_First_Biospecimen_ID.dna') %>%
    rbind(rna_ids %>%
            dplyr::rename('Kids_First_Biospecimen_ID' = 'Kids_First_Biospecimen_ID.rna')) %>%
    filter(Kids_First_Biospecimen_ID %in% bs_ids,
           sample_id %in% samples) %>%
    dplyr::select(Kids_First_Biospecimen_ID, sample_id) %>%
    unique()
}

# modify MAF
# determine which columns to keep 
keep_cols <- c("Chromosome",
               "Start_Position",
               "End_Position",
               "Strand",
               "Variant_Classification",
               "IMPACT",
               "Tumor_Sample_Barcode",
               "Hugo_Symbol",
               "HGVSp_Short",
               "Exon_Number",
               "CLIN_SIG")

# determine which small variant types to include
include_var_class <- c(
  "Missense_Mutation",
  "Frame_Shift_Del",
  "In_Frame_Ins",
  "Frame_Shift_Ins",
  "Splice_Site",
  "Nonsense_Mutation",
  "In_Frame_Del",
  "Nonstop_Mutation",
  "Translation_Start_Site"
)

maf_file <- maf_file %>%
  dplyr::filter(Variant_Classification %in% include_var_class,
                IMPACT %in% c('MODERATE', 'HIGH'),
                grepl('pathogenic|drug_response', CLIN_SIG)) %>%
  mutate(gene_symbol = Hugo_Symbol) %>%
  dplyr::select(Kids_First_Biospecimen_ID, gene_symbol) %>%
  dplyr::filter(gene_symbol %in% genes) %>%
  dplyr::filter(Kids_First_Biospecimen_ID %in% bs_ids) %>%
  unique()

# modify SV   
if(length(sv_file) > 0){
  sv_file <- sv_file %>%
    dplyr::filter(FILTER == 'PASS') %>%
    mutate(Kids_First_Biospecimen_ID = Kids.First.Biospecimen.ID.Tumor,
           gene_symbol = Gene.name) %>%
    dplyr::select(Kids_First_Biospecimen_ID, gene_symbol) %>%
    dplyr::filter(gene_symbol %in% genes) %>%
    dplyr::filter(Kids_First_Biospecimen_ID %in% bs_ids) %>%
    unique()
} else {
  sv_file <- data.frame(Kids_First_Biospecimen_ID = NA, gene_symbol = NA)
}

# modify CNV  
cnv_file <- cnv_file %>%
  mutate(gene_symbol = hgnc_symbol) %>%
  filter(gene_symbol %in% onco_genes & status %in% c("Gain", "Amplification") | 
           gene_symbol %in% tsg_genes & status %in% c("Loss", "Complete Loss")) %>%
  dplyr::select(Kids_First_Biospecimen_ID, gene_symbol) %>%
  dplyr::filter(gene_symbol %in% genes) %>%
  dplyr::filter(Kids_First_Biospecimen_ID %in% bs_ids) %>%
  unique()

# modify fusion  
if(length(fusion_file) > 0){
  fusion_file <- fusion_file %>%
    mutate(gene_symbol = strsplit(fusion_name, "_")) %>% 
    unnest(gene_symbol) %>% # split into new rows 
    dplyr::select(Kids_First_Biospecimen_ID, gene_symbol) %>%
    dplyr::filter(gene_symbol %in% genes) %>%
    dplyr::filter(Kids_First_Biospecimen_ID %in% bs_ids) %>%
    unique()
} else {
  fusion_file <- data.frame(Kids_First_Biospecimen_ID = NA, gene_symbol = NA)
}

# binary matrix of alterations
all_alterations <- rbind(maf_file, fusion_file, sv_file, cnv_file) %>%
  inner_join(histology_df, by = "Kids_First_Biospecimen_ID") %>%
  dplyr::select(-c(Kids_First_Biospecimen_ID)) %>%
  unique() 
all_alterations <- acast(all_alterations, gene_symbol~sample_id, fun.aggregate = length)
all_alterations <- as.data.frame(all_alterations)

# check if the sample of interest is present
# it may not have passed the filters
if(length(intersect(colnames(all_alterations), sample_id_interest)) == 0){
  all_alterations[,sample_id_interest] <- 0
}

# write file
write.table(all_alterations, file = mat_file, sep = "\t", col.names = NA, row.names = TRUE)
