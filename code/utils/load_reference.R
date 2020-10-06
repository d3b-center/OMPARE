#########################
# Load all reference data
#########################

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir('.git'))
source(file.path(root_dir, 'code', 'utils', 'define_directories.R'))

# GTEx Normals Brain TPM (1152 samples)
gtexData <- readRDS(file.path(ref_dir, 'GTEx', 'GTEx_Brain_TPM.RDS'))

# All PNOC008 patients (TPM matrix + clinical)
pnoc008.data <- readRDS(file.path(ref_dir, 'PNOC008', 'PNOC008_TPM_matrix.RDS'))
pnoc008.clinData <- readRDS(file.path(ref_dir, 'PNOC008', 'PNOC008_clinData.RDS'))

# PBTA specific mRNA data (TPM)
pbta.full <- readRDS(file.path(ref_dir, 'PBTA', 'pbta-gene-expression-rsem-tpm-collapsed.polya.stranded.rds'))
pbta.mat <- pbta.full
pbta.clinData <- read.delim(file.path(ref_dir, 'PBTA', 'pbta-histologies.tsv'), stringsAsFactors = F)
pbta.survData <- pbta.clinData %>%
  filter(experimental_strategy == 'RNA-Seq') %>%
  mutate(sample_barcode = Kids_First_Biospecimen_ID) %>%
  dplyr::select(sample_barcode, sample_id, OS_days, OS_status) %>%
  filter(!is.na(OS_status)) %>%
  mutate(OS_status = ifelse(OS_status == 'DECEASED', 1, 0))

# PBTA specific CNV data
pbta.cnv <- data.table::fread(file.path(ref_dir, 'PBTA', 'pbta-cnv-controlfreec.tsv.gz'))

# TCGA GBM specific data (TPM)
tcga.gbm.mat <- readRDS(file.path(ref_dir, 'TCGA', 'TCGA_GBM_matrix_TPM.RDS'))
tcga.gbm.clinData <- readRDS(file.path(ref_dir, 'TCGA', 'TCGA_GBM_clinData.RDS'))
tcga.gbm.survData <- tcga.gbm.clinData %>%
  filter(overall_survival_time_in_days != 'unavailable') %>%
  mutate(OS_status = as.numeric(vital_status),
         OS_days = as.numeric(overall_survival_time_in_days)) %>%
  dplyr::select(sample_barcode, OS_days, OS_status)

# Cancer Genes
cancerGenes <- readRDS(file.path(ref_dir, 'cancer_gene_list.rds'))
tsgGenes <- read.delim(file.path(ref_dir, 'Human_TSGs.txt'), stringsAsFactors = F)

# CNV analysis (chr coordinates to gene symbol map)
chrMap <- read.delim(file.path(ref_dir, 'mart_export_genechr_mapping.txt'), stringsAsFactors = F, check.names = F)
colnames(chrMap) <- c("hgnc_symbol", "gene_start", "gene_end", "chromosome")

# Network analysis
geneMania <- read.delim(file.path(ref_dir, 'GeneManiaNetwork.txt'), stringsAsFactors =F)

# HGG-specific genes
diseaseSpecificFields <- read.delim(file.path(ref_dir, 'DiseaseSpecificFields.txt'))

# Gene sets (c2 reactome)
geneSet <- getGmt(file.path(ref_dir, 'mSigDB', 'c2.cp.reactome.v6.0.symbols.gmt'), collectionType = BroadCollection(), geneIdType = SymbolIdentifier())
geneSet <- geneIds(geneSet)
geneSetTS <- stack(geneSet)

# Mutational signatures
signatures <- readAlexandrovSignatures(file.path(ref_dir, 'signatures_probabilities.txt'))

# Germline markers
germlineMarkers <- readRDS(file.path(ref_dir, 'germline_markers_list.rds'))

# TMB from PBTA and TCGA
pedTMB <- data.table::fread(file.path(ref_dir, 'pbta-TMBscores_withdiseastype.txt'))
adultTMB <- data.table::fread(file.path(ref_dir, 'TCGA_diseasetypes_and_samples_TMBscores.txt'))
TMBFileBED <- data.table::fread(file.path(ref_dir, 'xgen-exome-research-panel-targets_hg38_ucsc_liftover.100bp_padded.sort.merged.bed'))
colnames(TMBFileBED)  <- c('chr', 'start', 'end')

# genelist for heatmaps
genelist.heatmap <- read.delim(file.path(ref_dir, '2020-03-30_Glioma_GeneList.txt'), stringsAsFactors = F)

