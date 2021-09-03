#########################
# Load all reference data
#########################

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir('.git'))
data_dir <- file.path(root_dir, "data")

# gencode reference
gencode_v27 <- read.delim(file.path(data_dir, 'pnoc008', 'gencode.v27.primary_assembly.annotation.txt'))
gencode_v27_pc <- gencode_v27 %>%
  filter(biotype == "protein_coding")

# GTEx Normals Brain TPM (1152 samples)
gtex_brain_tpm <- readRDS(file.path(data_dir, 'gtex', 'gtex_brain_tpm.rds'))
gtex_brain_counts <- readRDS(file.path(data_dir, 'gtex', 'gtex_brain_counts.rds'))

# All PNOC008 patients (TPM matrix + clinical)
pnoc008_tpm <- readRDS(file.path(data_dir, 'pnoc008', 'pnoc008_tpm_matrix.rds'))
pnoc008_clinical <- readRDS(file.path(data_dir, 'pnoc008', 'pnoc008_clinical.rds'))

# PBTA specific mRNA data (TPM)
pbta_full_tpm <- readRDS(file.path(data_dir, 'pbta', 'pbta-gene-expression-rsem-tpm-collapsed.polya.stranded.rds'))
pbta_tpm <- pbta_full_tpm
pbta_clinical <- read.delim(file.path(data_dir, 'pbta', 'pbta-histologies.tsv'), stringsAsFactors = F)

# survival needs to be restricted to HGAT samples only
pbta_survival <- pbta_clinical %>%
  filter(experimental_strategy == 'RNA-Seq',
         short_histology == "HGAT",
         !is.na(OS_status)) %>%
  mutate(subject_id = Kids_First_Biospecimen_ID,
         OS_status = ifelse(OS_status == 'DECEASED', 1, 0)) %>%
  dplyr::select(subject_id, OS_days, OS_status) 

# TCGA GBM specific data (TPM)
tcga_gbm_tpm <- readRDS(file.path(data_dir, 'tcga', 'tcga_gbm_tpm_matrix.rds'))
tcga_gbm_clinical <- readRDS(file.path(data_dir, 'tcga', 'tcga_gbm_clinical.rds'))
tcga_gbm_survival <- tcga_gbm_clinical %>%
  filter(overall_survival_time_in_days != 'unavailable') %>%
  mutate(subject_id = sample_barcode,
         OS_status = as.numeric(vital_status),
         OS_days = as.numeric(overall_survival_time_in_days)) %>%
  dplyr::select(subject_id, OS_days, OS_status)

# Cancer Genes (annoFuse)
cancer_genes <- readRDS(file.path(data_dir, 'cancer_gene_list.rds'))

# CNV analysis (chr coordinates to gene symbol map)
chr_map <- read.delim(file.path(data_dir, 'mart_export_genechr_mapping.txt'), stringsAsFactors = F, check.names = F)
colnames(chr_map) <- c("hgnc_symbol", "gene_start", "gene_end", "chromosome")

# Network analysis
# gene_mania <- read.delim(file.path(data_dir, 'GeneManiaNetwork.txt'), stringsAsFactors =F)

# Gene sets (c2 reactome)
gene_set <- getGmt(file.path(data_dir, 'msigdb', 'c2.cp.reactome.v6.0.symbols.gmt'), collectionType = BroadCollection(), geneIdType = SymbolIdentifier())
gene_set <- geneIds(gene_set)
gene_set_ts <- stack(gene_set)

# genelist for heatmaps
# genelist_heatmap <- read.delim(file.path(data_dir, '2020-03-30_Glioma_GeneList.txt'), stringsAsFactors = F)

# cancer hotspot (https://www.cancerhotspots.org/#/download)
cancer_hotspots_v2  <- readRDS(file.path(data_dir, 'cancer_hotspots_v2.rds'))
