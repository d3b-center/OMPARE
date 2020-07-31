#########################
# Load all reference data
#########################

# GTEx Normals TPM (7863 samples)
gtexData <- readRDS("data/Reference/GTEx/gtex_normals_TPM.RDS")

# All PNOC008 patients (TPM matrix + clinical)
pnoc008.data <- readRDS('data/Reference/PNOC008/PNOC008_TPM_matrix.RDS')
pnoc008.clinData <- readRDS('data/Reference/PNOC008/PNOC008_clinData.RDS')

# PBTA specific mRNA data (TPM)
pbta.full <- readRDS('data/Reference/PBTA/pbta-gene-expression-rsem-tpm-collapsed.polya.stranded.rds')
pbta.mat <- pbta.full
pbta.clinData <- read.delim("data/Reference/PBTA/pbta-histologies.tsv", stringsAsFactors = F)
pbta.survData <- pbta.clinData %>%
  filter(experimental_strategy == "RNA-Seq") %>%
  mutate(sample_barcode = Kids_First_Biospecimen_ID) %>%
  dplyr::select(sample_barcode, sample_id, OS_days, OS_status) %>%
  filter(!is.na(OS_status)) %>%
  mutate(OS_status = ifelse(OS_status == "DECEASED", 1, 0))

# PBTA specific CNV data
pbta.cnv <- data.table::fread('data/Reference/PBTA/pbta-cnv-controlfreec.tsv.gz')

# TCGA GBM specific data (TPM)
tcga.gbm.mat <- readRDS('data/Reference/TCGA/TCGA_GBM_matrix_TPM.RDS')
tcga.gbm.clinData <- readRDS('data/Reference/TCGA/TCGA_GBM_clinData.RDS')
tcga.gbm.survData <- tcga.gbm.clinData %>%
  filter(overall_survival_time_in_days != "unavailable") %>%
  mutate(OS_status = as.numeric(vital_status),
         OS_days = as.numeric(overall_survival_time_in_days)) %>%
  dplyr::select(sample_barcode, OS_days, OS_status)

# Cancer Genes
cancerGenes <- read.delim("data/Reference/CancerGeneList.tsv", stringsAsFactors = F)
cancerGenes <- cancerGenes %>%
  filter(Gene_Symbol != "") %>%
  dplyr::select(-Count) %>%
  gather(key = "file", value = "type", -Gene_Symbol) %>%
  mutate(type = file)
geneListRef <- read.delim("data/Reference/genelistreference.txt", stringsAsFactors = F)
geneListRef <- subset(geneListRef, type == "TumorSuppressorGene" | type == "CosmicCensus" | type == "Oncogene")
cancerGenes <- rbind(cancerGenes, geneListRef)
rm(geneListRef)
tsgGenes <- read.delim("data/Reference/Human_TSGs.txt", stringsAsFactors = F)

# CNV analysis (chr coordinates to gene symbol map)
chrMap <- read.delim("data/Reference/mart_export_genechr_mapping.txt", stringsAsFactors =F)

# Network analysis
geneMania <- read.delim("data/Reference/GeneManiaNetwork.txt", stringsAsFactors =F)

# HGG-specific genes
diseaseSpecificFields <- read.delim("data/Reference/DiseaseSpecificFields.txt")

# Gene sets (KEGG)
geneSet <- getGmt('data/Reference/mSigDB/c2.cp.kegg.v7.1.symbols.gmt', collectionType = BroadCollection(), geneIdType = SymbolIdentifier())
geneSet <- geneIds(geneSet)
geneSetTS <- stack(geneSet)

# Mutational signatures
signatures <- readAlexandrovSignatures("data/Reference/signatures_probabilities.txt")

# Drug info
dgidb <- read.delim("data/Reference/DGIdb.txt", stringsAsFactors = F)

# Germline markers
pharmacogenomics.genes <- read.delim("data/Reference/Pharmacogenomics_Genes.list", stringsAsFactors = F, header = F)
pharmacogenomics.genes <- data.frame("Gene" = pharmacogenomics.genes$V1, "Class" = "Pharmacogenomics")
chop.panel.genes <- read.delim("data/Reference/CHOP_Additional_Cancer_Genes.list", stringsAsFactors = F, header = F)
chop.panel.genes <-  data.frame("Gene" = chop.panel.genes$V1, "Class" = "CHOP Panel")
acmg.genes <- read.delim("data/Reference/ACMG_22_Cancer_Genes.list", stringsAsFactors = F, header = F)
acmg.genes <-  data.frame("Gene" = acmg.genes$V1, "Class" = "ACMG")
germlineMarkers <- rbind(pharmacogenomics.genes,  acmg.genes, chop.panel.genes)

# TMB from PBTA and TCGA
pedTMB <- data.table::fread("data/Reference/pbta-TMBscores_withdiseastype.txt")
adultTMB <- data.table::fread("data/Reference/TCGA_diseasetypes_and_samples_TMBscores.txt")
TMBFileBED <- data.table::fread("data/Reference/xgen-exome-research-panel-targets_hg38_ucsc_liftover.100bp_padded.sort.merged.bed")
colnames(TMBFileBED)  <- c("chr", "start", "end")

# genelist for heatmaps
genelist.heatmap <- read.delim('data/Reference/2020-03-30_Glioma_GeneList.txt', stringsAsFactors = F)

