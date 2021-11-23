# Author: Run Jin
#
# GSNCA analysis comparing upper and lower quantile of gene expressions in each disease
# BiocManager::install("GSAR")
# BiocManager::install("GSVAdata")
# BiocManager::install("DGCA")
# BiocManager::install("EGSEA")
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("GSAR"))
suppressPackageStartupMessages(library("org.Hs.eg.db"))
suppressPackageStartupMessages(library("EGSEA"))
suppressPackageStartupMessages(library("DGCA"))

# arguments
option_list <- list(
  make_option(c("--patient"), type = "character",
              help = "Patient identifier (PNOC008-22, C3342894...)")
)
opt <- parse_args(OptionParser(option_list = option_list))
patient <- opt$patient

#### Define Directories --------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
module_dir <- file.path(root_dir, "code", "gsnca_analysis")
patient_dir <- file.path(root_dir, "results", patient)

input_dir <- file.path(patient_dir, "output", "transcriptomically_similar_analysis")
output_dir <- file.path(patient_dir, "output", "gsnca_analysis")
dir.create(output_dir, showWarnings = F, recursive = T)

#### Source utils functions --------------------------------------------------------
source(file.path(module_dir, "utils", "gsnca_calc.R"))

#### Read in files necessary for analyses --------------------------------------
# UMAP results
transciptomically_similar_pediatric <- readRDS(file.path(input_dir, "transciptomically_similar_pediatric.rds"))
transciptomically_similar_adult <- readRDS(file.path(input_dir, "transciptomically_similar_adult.rds"))

# gencode reference
gencode_v27 <- read.delim(file.path(data_dir, 'pnoc008', 'gencode.v27.primary_assembly.annotation.txt'))
gencode_v27_pc <- gencode_v27 %>%
  filter(biotype == "protein_coding")

# pnoc008 samples 
pnoc008_histology <- readRDS(file.path(data_dir, 'pnoc008', 'pnoc008_clinical.rds'))
pnoc008_tpm<- readRDS(file.path(data_dir, 'pnoc008', 'pnoc008_tpm_matrix.rds')) %>% 
  tibble::rownames_to_column("geneID")

# GTex Brain TPM
gtex_brain_clinical <- readRDS(file.path(data_dir, 'gtex', 'gtex_brain_clinical.rds'))
gtex_brain_tpm <- readRDS(file.path(data_dir, "gtex", "gtex_brain_tpm.rds"))
gtex_brain_tpm <- gtex_brain_tpm[grep("^HIST", rownames(gtex_brain_tpm), invert = T),]
gtex_brain_tpm <- gtex_brain_tpm[rownames(gtex_brain_tpm) %in% gencode_v27_pc$gene_symbol,]

###### PBTA
# clinical
pbta_clinical <- read.delim(file.path(data_dir, 'pbta', 'pbta-histologies.tsv'), stringsAsFactors = F)

# expression  (polyA + stranded combined data collapsed to gene symbols)
pbta_full_tpm <- readRDS(file.path(data_dir, 'pbta','pbta-gene-expression-rsem-tpm-collapsed.polya.stranded.rds'))
pbta_full_tpm <- pbta_full_tpm[grep("^HIST", rownames(pbta_full_tpm), invert = T),]
pbta_full_tpm <- pbta_full_tpm[rownames(pbta_full_tpm) %in% gencode_v27_pc$gene_symbol,]

# filter pbta clinical to HGAT
pbta_clinical_hgg <- pbta_clinical %>%
  filter(experimental_strategy == "RNA-Seq",
         short_histology == "HGAT",
         Kids_First_Biospecimen_ID %in% colnames(pbta_full_tpm))

pbta_hgg_tpm <- pbta_full_tpm[,colnames(pbta_full_tpm) %in% pbta_clinical_hgg$Kids_First_Biospecimen_ID]

#### Get expression matrix of POI + 20 similar samples --------------------------------------
pnoc008_pbta_tpm <- pbta_full_tpm %>% 
  tibble::rownames_to_column("geneID") %>%
  full_join(pnoc008_tpm) %>%
  replace(is.na(.), 0)

similar_subject_ids<- transciptomically_similar_pediatric %>% 
  pull(subject_id)
stopifnot(all_of(similar_subject_ids) %in% colnames(pnoc008_pbta_tpm))
  
pnoc008_similar_tpm <- pnoc008_pbta_tpm %>%
  dplyr::select(geneID, all_of(patient), all_of(similar_subject_ids)) %>%
  tibble::column_to_rownames("geneID")

#### Prepare background/cohort comparison expression matrix -------------------------------
# gtex matrix 
gtex_brain_tpm_filtered <- filter_low_expr_df(gtex_brain_tpm)

# the rest of HGG
pbta_hgg_tpm_rest_filtered <- pbta_hgg_tpm %>%
  dplyr::select(-all_of(similar_subject_ids[similar_subject_ids %in% colnames(pbta_hgg_tpm)])) %>%
  filter_low_expr_df()

# the rest of PBTA 
pbta_full_tpm_rest_filtered <- pbta_full_tpm %>%
  dplyr::select(-all_of(similar_subject_ids[similar_subject_ids %in% colnames(pbta_full_tpm)])) %>%
  filter_low_expr_df()

#### Run GSNCA and output results in text files ----------------------
gsnca_analysis_plot(pnoc008_similar_tpm, gtex_brain_tpm_filtered, "GTEx", top_bar=20, top_net=5)

gsnca_analysis_plot(pnoc008_similar_tpm, pbta_hgg_tpm_rest_filtered, "PBTA_HGG", top_bar=20, top_net=5)

gsnca_analysis_plot(pnoc008_similar_tpm, pbta_full_tpm_rest_filtered, "PBTA_all", top_bar=20, top_net=5)


