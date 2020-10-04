# Author: Komal S. Rathi
# Date: 04/25/2020
# Function: Up/Down pathways for each PBTA sample, compare to rest of PBTA (1) and GTEx (2) and TCGA GBM (3) and PNOC008 (4)
# do this once and read in for tabulate pathways (Page 8 of report)

# Function to return all results from RNA-Seq Analysis
suppressPackageStartupMessages(library(GSEABase))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(xlsx))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(doMC))
registerDoMC(cores = 4)

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source function for RNA-seq diffexpr & pathway analysis
source(file.path(utils_dir, "rnaseq_analysis_accessory.R"))

# Dataset1: GTex Brain
gtexBrain <- readRDS(file.path(ref_dir, 'GTEx', 'GTEx_Brain_clinical.RDS'))
gtexData <- readRDS(file.path(ref_dir, "GTEx", "GTEx_Brain_TPM.RDS"))
gtexData <- gtexData[grep("^HIST", rownames(gtexData), invert = T),]

# Dataset2: TCGA GBM
tcgaGBMclin <- readRDS(file.path(ref_dir, 'TCGA', 'TCGA_GBM_clinData.RDS'))
tcgaGBMData <- readRDS(file.path(ref_dir, 'TCGA', 'TCGA_GBM_matrix_TPM.RDS'))
tcgaGBMData <- tcgaGBMData[grep("^HIST", rownames(tcgaGBMData), invert = T),]

# Dataset3: PBTA (polyA + corrected stranded n = 1028)
# clinical
pbta.hist <- read.delim(file.path(ref_dir, 'PBTA', 'pbta-histologies.tsv'), stringsAsFactors = F)
pbta.hist <- pbta.hist %>%
  filter(experimental_strategy == "RNA-Seq",
         short_histology == "HGAT")

# expression  (polyA + stranded combined TPM data collapsed to gene symbols)
pbta.full <- readRDS(file.path(ref_dir, 'PBTA','pbta-gene-expression-rsem-tpm-collapsed.polya.stranded.rds'))
pbta.full <- pbta.full[grep("^HIST", rownames(pbta.full), invert = T),]

# Dataset4: PBTA (polyA + corrected stranded HGG n = 186)
pbta.hgg <- pbta.full[,colnames(pbta.full) %in% pbta.hist$Kids_First_Biospecimen_ID]

# Dataset5: PNOC008
pnoc008 <- readRDS(file.path(ref_dir, 'PNOC008', 'PNOC008_TPM_matrix.RDS'))
pnoc008 <- pnoc008[grep("^HIST", rownames(pnoc008), invert = T),]

# Cancer Genes
cancerGenes <- readRDS(file.path(ref_dir, 'cancer_gene_list.rds'))

# Genesets (c2 reactome)
geneSet <- getGmt(file.path(ref_dir, 'mSigDB', 'c2.cp.reactome.v6.0.symbols.gmt'), collectionType = BroadCollection(), geneIdType = SymbolIdentifier())
geneSet <- geneIds(geneSet)

# input data
res.pbta <- melt(as.matrix(pbta.full), value.name = "TPM", varnames = c("Gene", "Sample"))
res.tcga <- melt(as.matrix(tcgaGBMData), value.name = "TPM", varnames = c("Gene", "Sample"))
res.pnoc008 <- melt(as.matrix(pnoc008), value.name = "TPM", varnames = c("Gene", "Sample"))

# create output directory
gsea.dir <- file.path(ref_dir, 'GSEA')
dir.create(gsea.dir, showWarnings = F, recursive = T)

# overwrite PNOC008 comparisons 
# PNOC008 vs GTEx Brain
GTExBrain <-  plyr::dlply(res.pnoc008, .variables = "Sample", .fun = function(x) runRNASeqAnalysis(exp.data = x, refData = gtexData, comparison = paste0("GTExBrain_", ncol(gtexData))), .parallel = TRUE)
saveRDS(GTExBrain, file = file.path(ref_dir, 'GSEA', 'PNOC008_vs_GTExBrain.RDS'))

# PNOC008 vs PBTA
PBTA_All <- plyr::dlply(res.pnoc008, .variables = "Sample", .fun = function(x) runRNASeqAnalysis(exp.data = x, refData = pbta.full, comparison = paste0("PBTA_All_", ncol(pbta.full))), .parallel = TRUE)
saveRDS(PBTA_All, file = file.path(ref_dir, 'GSEA', 'PNOC008_vs_PBTA.RDS'))

# PNOC008 vs PBTA HGG
PBTA_HGG <- plyr::dlply(res.pnoc008, .variables = "Sample", .fun = function(x) runRNASeqAnalysis(exp.data = x, refData = pbta.hgg, comparison = paste0("PBTA_HGG_", ncol(pbta.hgg))), .parallel = TRUE)
saveRDS(PBTA_HGG, file = file.path(ref_dir, 'GSEA', 'PNOC008_vs_PBTA_HGG.RDS'))

# PNOC008 vs TCGA GBM
TCGA_GBM <- plyr::dlply(res.pnoc008, .variables = "Sample", .fun = function(x) runRNASeqAnalysis(exp.data = x, refData = tcgaGBMData, comparison = paste0("TCGA_GBM_", ncol(tcgaGBMData))), .parallel = TRUE)
saveRDS(TCGA_GBM, file = file.path(ref_dir, 'GSEA', 'PNOC008_vs_TCGA_GBM.RDS'))

# PBTA comparisons only need to be run once
# PBTA vs GTEx Brain
if(!file.exists(file.path(ref_dir, 'GSEA', 'PBTA_vs_GTExBrain.RDS'))){
  GTExBrain <-  plyr::dlply(res.pbta, .variables = "Sample", .fun = function(x) runRNASeqAnalysis(exp.data = x, refData = gtexData, comparison = paste0("GTExBrain_", ncol(gtexData))), .parallel = TRUE)
  saveRDS(GTExBrain, file = file.path(ref_dir, 'GSEA', 'PBTA_vs_GTExBrain.RDS'))
}

# PBTA vs PBTA
if(!file.exists(file.path(ref_dir, 'GSEA', 'PBTA_vs_PBTA.RDS'))){
  PBTA_All <- plyr::dlply(res.pbta, .variables = "Sample", .fun = function(x) runRNASeqAnalysis(exp.data = x, refData = pbta.full, comparison = paste0("PBTA_All_", ncol(pbta.full))), .parallel = TRUE)
  saveRDS(PBTA_All, file = file.path(ref_dir, 'GSEA', 'PBTA_vs_PBTA.RDS'))
}

# PBTA vs PBTA HGG
if(!file.exists(file.path(ref_dir, 'GSEA', 'PBTA_vs_PBTAHGG.RDS'))){
  PBTA_HGG <- plyr::dlply(res.pbta, .variables = "Sample", .fun = function(x) runRNASeqAnalysis(exp.data = x, refData = pbta.hgg, comparison = paste0("PBTA_HGG_", ncol(pbta.hgg))), .parallel = TRUE)
  saveRDS(PBTA_HGG, file = file.path(ref_dir, 'GSEA', 'PBTA_vs_PBTAHGG.RDS'))
}

# TCGA comparisons only need to be run once
# TCGA GBM vs GTEx Brain
if(!file.exists(file.path(ref_dir, 'GSEA', 'TCGA_GBM_vs_GTExBrain.RDS'))){
  TCGA.vs.GTExBrain <-  plyr::dlply(res.tcga, .variables = "Sample", .fun = function(x) runRNASeqAnalysis(exp.data = x, refData = gtexData, comparison = paste0("GTExBrain_", ncol(gtexData))), .parallel = TRUE)
  saveRDS(TCGA.vs.GTExBrain, file = file.path(ref_dir, 'GSEA', 'TCGA_GBM_vs_GTExBrain.RDS'))
}

# TCGA GBM vs TCGA GBM
if(!file.exists(file.path(ref_dir, 'GSEA', 'TCGA_GBM_vs_TCGA_GBM.RDS'))){
  TCGA.vs.TCGA <-  plyr::dlply(res.tcga, .variables = "Sample", .fun = function(x) runRNASeqAnalysis(exp.data = x, refData = tcgaGBMData, comparison = paste0("TCGA_GBM_", ncol(tcgaGBMData))), .parallel = TRUE)
  saveRDS(TCGA.vs.TCGA, file = file.path(ref_dir, 'GSEA', 'TCGA_GBM_vs_TCGA_GBM.RDS'))
}
