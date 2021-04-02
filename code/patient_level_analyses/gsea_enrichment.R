# Author: Komal S. Rathi
# Date: 04/25/2020
# Function: Up/Down pathways for each PBTA sample, compare to rest of PBTA (1) and GTEx (2) and TCGA GBM (3) and PNOC008 (4)
# do this once and read in for tabulate pathways (Page 8 of report)

# Function to return all results from RNA-Seq Analysis
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(GSEABase))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(doMC))
suppressPackageStartupMessages(library(optparse))
registerDoMC(cores = 2)

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

option_list <- list(
  make_option(c("-p", "--patient"), type = "character",
              help = "Patient identifier for e.g. PNOC008-1, PNOC008-10 etc")
)
opt <- parse_args(OptionParser(option_list = option_list))
sample_to_add <- opt$patient

# source function for RNA-seq diffexpr & pathway analysis
source(file.path(patient_level_analyses_utils, "rnaseq_analysis_edgeR.R"))

# gencode reference
gencode_v27 <- read.delim(file.path(ref_dir, 'pnoc008', 'gencode.v27.primary_assembly.annotation.txt'))
gencode_v27_pc <- gencode_v27 %>%
  filter(biotype == "protein_coding")

# Dataset1: GTex Brain
gtex_brain_clinical <- readRDS(file.path(ref_dir, 'gtex', 'gtex_brain_clinical.rds'))
gtex_brain_tpm <- readRDS(file.path(ref_dir, "gtex", "gtex_brain_tpm.rds"))
gtex_brain_tpm <- gtex_brain_tpm[grep("^HIST", rownames(gtex_brain_tpm), invert = T),]
gtex_brain_tpm <- gtex_brain_tpm[rownames(gtex_brain_tpm) %in% gencode_v27_pc$gene_symbol,]

gtex_brain_counts <- readRDS(file.path(ref_dir, "gtex", "gtex_brain_counts.rds"))
gtex_brain_counts <- gtex_brain_counts[grep("^HIST", rownames(gtex_brain_counts), invert = T),]
gtex_brain_counts <- gtex_brain_counts[rownames(gtex_brain_counts) %in% gencode_v27_pc$gene_symbol,]

# Dataset2: TCGA GBM
tcga_gbm_clinical <- readRDS(file.path(ref_dir, 'tcga', 'tcga_gbm_clinical.rds'))
tcga_gbm_tpm <- readRDS(file.path(ref_dir, 'tcga', 'tcga_gbm_tpm_matrix.rds'))
tcga_gbm_tpm <- tcga_gbm_tpm[grep("^HIST", rownames(tcga_gbm_tpm), invert = T),]
tcga_gbm_tpm <- tcga_gbm_tpm[rownames(tcga_gbm_tpm) %in% gencode_v27_pc$gene_symbol,]

tcga_gbm_counts <- readRDS(file.path(ref_dir, 'tcga', 'tcga_gbm_counts_matrix.rds'))
tcga_gbm_counts <- tcga_gbm_counts[grep("^HIST", rownames(tcga_gbm_counts), invert = T),]
tcga_gbm_counts <- tcga_gbm_counts[rownames(tcga_gbm_counts) %in% gencode_v27_pc$gene_symbol,]

# Dataset3: PBTA (polyA + corrected stranded n = 1028)
# clinical
pbta_clinical <- read.delim(file.path(ref_dir, 'pbta', 'pbta-histologies.tsv'), stringsAsFactors = F)
pbta_clinical <- pbta_clinical %>%
  filter(experimental_strategy == "RNA-Seq",
         short_histology == "HGAT")

# expression  (polyA + stranded combined data collapsed to gene symbols)
pbta_full_tpm <- readRDS(file.path(ref_dir, 'pbta','pbta-gene-expression-rsem-tpm-collapsed.polya.stranded.rds'))
pbta_full_tpm <- pbta_full_tpm[grep("^HIST", rownames(pbta_full_tpm), invert = T),]
pbta_full_tpm <- pbta_full_tpm[rownames(pbta_full_tpm) %in% gencode_v27_pc$gene_symbol,]

pbta_full_counts <- readRDS(file.path(ref_dir, 'pbta','pbta-gene-expression-rsem-counts-collapsed.polya.stranded.rds'))
pbta_full_counts <- pbta_full_counts[grep("^HIST", rownames(pbta_full_counts), invert = T),]
pbta_full_counts <- pbta_full_counts[rownames(pbta_full_counts) %in% gencode_v27_pc$gene_symbol,]

# Dataset4: PBTA (polyA + corrected stranded HGG n = 189)
pbta_hgg_tpm <- pbta_full_tpm[,colnames(pbta_full_tpm) %in% pbta_clinical$Kids_First_Biospecimen_ID]
pbta_hgg_counts <- pbta_full_counts[,colnames(pbta_full_counts) %in% pbta_clinical$Kids_First_Biospecimen_ID]

# Dataset5: PNOC008
pnoc008_tpm <- readRDS(file.path(ref_dir, 'pnoc008', 'pnoc008_tpm_matrix.rds'))
pnoc008_tpm <- pnoc008_tpm[grep("^HIST", rownames(pnoc008_tpm), invert = T),]
pnoc008_tpm <- pnoc008_tpm[rownames(pnoc008_tpm) %in% gencode_v27_pc$gene_symbol,]

pnoc008_counts <- readRDS(file.path(ref_dir, 'pnoc008', 'pnoc008_counts_matrix.rds'))
pnoc008_counts <- pnoc008_counts[grep("^HIST", rownames(pnoc008_counts), invert = T),]
pnoc008_counts <- pnoc008_counts[rownames(pnoc008_counts) %in% gencode_v27_pc$gene_symbol,]

# Cancer Genes
cancer_genes <- readRDS(file.path(ref_dir, 'cancer_gene_list.rds'))

# Genesets (c2 reactome)
gene_set <- getGmt(file.path(ref_dir, 'msigdb', 'c2.cp.reactome.v6.0.symbols.gmt'), collectionType = BroadCollection(), geneIdType = SymbolIdentifier())
gene_set <- geneIds(gene_set)

# DSigDb geneset
dsigdb_geneset <- getGmt(file.path(ref_dir, 'dsigdb', 'DSigDB_All.gmt'))
dsigdb_geneset <- geneIds(dsigdb_geneset)
dsigdb_geneset <- dsigdb_geneset[grep('_TTD_|_BOSS$|_CTD_', names(dsigdb_geneset), value = T, invert = T)]

# input data
pbta_full_tpm_melt <- melt(as.matrix(pbta_full_tpm), value.name = "TPM", varnames = c("Gene", "Sample"))
tcga_gbm_tpm_melt <- melt(as.matrix(tcga_gbm_tpm), value.name = "TPM", varnames = c("Gene", "Sample"))
pnoc008_tpm_melt <- melt(as.matrix(pnoc008_tpm), value.name = "TPM", varnames = c("Gene", "Sample"))
pbta_full_counts_melt <- melt(as.matrix(pbta_full_counts), value.name = "Counts", varnames = c("Gene", "Sample"))
tcga_gbm_counts_melt <- melt(as.matrix(tcga_gbm_counts), value.name = "Counts", varnames = c("Gene", "Sample"))
pnoc008_counts_melt <- melt(as.matrix(pnoc008_counts), value.name = "Counts", varnames = c("Gene", "Sample"))

# create output directory
gsea.dir <- file.path(ref_dir, 'gsea')
pnoc008.dir <- file.path(ref_dir, 'pnoc008')
dir.create(gsea.dir, showWarnings = F, recursive = T)
dir.create(pnoc008.dir, showWarnings = F, recursive = T)

# function to read existing gsea output, process and add current subject to it if not already present
gsea_enrichment <- function(sample_to_add, gsea_output, exp.data.counts, exp.data.tpm, refData.counts, gene_set, comparison){
  x <- exp.data.counts %>%
    filter(Sample %in% sample_to_add)
  if(file.exists(gsea_output)){
    existing_gsea_output <- readRDS(gsea_output)
  } else {
    existing_gsea_output <- list()
  }
  
  # run only if sample does not exist
  if(length(intersect(names(existing_gsea_output), sample_to_add)) == 0){
    res <- run_rnaseq_analysis_edger(exp.data.counts = x, 
                                     exp.data.tpm = exp.data.tpm, 
                                     refData.counts = refData.counts, 
                                     gene_set = gene_set, 
                                     comparison = comparison)
    existing_gsea_output[[sample_to_add]] <- res
    saveRDS(existing_gsea_output, file = gsea_output)
  }
}

# pnoc008 comparisons need to be updated with each new patient
# pnoc008 vs gtex brain
gsea_output <- file.path(ref_dir, 'gsea', 'pnoc008_vs_gtex_brain.rds')
gsea_enrichment(sample_to_add = sample_to_add, 
                gsea_output = gsea_output, 
                exp.data.counts = pnoc008_counts_melt, 
                exp.data.tpm = pnoc008_tpm_melt,
                refData.counts = gtex_brain_counts, 
                gene_set = gene_set, 
                comparison = paste0("GTExBrain_", ncol(gtex_brain_counts)))

# function to merge degs for pnoc008 vs gtex brain
merge_deg_list <- function(x){
  degs <- x$genes
  degs <- degs %>%
    dplyr::select(genes, diff_expr)
  return(degs)
}
pnoc008_vs_gtex_brain <- readRDS(gsea_output)
pnoc008_deg <- sapply(pnoc008_vs_gtex_brain, FUN = function(x) merge_deg_list(x = x), simplify = F, USE.NAMES = T)
pnoc008_deg <- data.table::rbindlist(pnoc008_deg, idcol = 'sample_name', use.names = T)
saveRDS(pnoc008_deg, file = file.path(pnoc008.dir, "pnoc008_vs_gtex_brain_degs.rds"))

# pnoc008 vs gtex brain (dsigdb)
gsea_output <- file.path(ref_dir, 'dsigdb', 'pnoc008_vs_gtex_brain.rds')
gsea_enrichment(sample_to_add = sample_to_add, 
                gsea_output = gsea_output, 
                exp.data.counts = pnoc008_counts_melt, 
                exp.data.tpm = pnoc008_tpm_melt,
                refData.counts = gtex_brain_counts, 
                gene_set = dsigdb_geneset, 
                comparison = paste0("GTExBrain_", ncol(gtex_brain_counts)))

# pnoc008 vs pbta
gsea_output <- file.path(ref_dir, 'gsea', 'pnoc008_vs_pbta.rds')
gsea_enrichment(sample_to_add = sample_to_add, 
                gsea_output = gsea_output, 
                exp.data.counts = pnoc008_counts_melt, 
                exp.data.tpm = pnoc008_tpm_melt,
                refData.counts = pbta_full_counts, 
                gene_set = gene_set, 
                comparison = paste0("PBTA_All_", ncol(pbta_full_counts)))


# pnoc008 vs pbta (dsigdb)
gsea_output <- file.path(ref_dir, 'dsigdb', 'pnoc008_vs_pbta.rds')
gsea_enrichment(sample_to_add = sample_to_add, 
                gsea_output = gsea_output, 
                exp.data.counts = pnoc008_counts_melt, 
                exp.data.tpm = pnoc008_tpm_melt,
                refData.counts = pbta_full_counts, 
                gene_set = dsigdb_geneset, 
                comparison = paste0("PBTA_All_", ncol(pbta_full_counts)))


# pnoc008 vs pbta hgg
gsea_output <- file.path(ref_dir, 'gsea', 'pnoc008_vs_pbta_hgg.rds')
gsea_enrichment(sample_to_add = sample_to_add, 
                gsea_output = gsea_output, 
                exp.data.counts = pnoc008_counts_melt, 
                exp.data.tpm = pnoc008_tpm_melt,
                refData.counts = pbta_hgg_counts, 
                gene_set = gene_set, 
                comparison = paste0("PBTA_HGG_", ncol(pbta_hgg_counts)))


# pnoc008 vs pbta hgg (dsigdb)
gsea_output <- file.path(ref_dir, 'dsigdb', 'pnoc008_vs_pbta_hgg.rds')
gsea_enrichment(sample_to_add = sample_to_add, 
                gsea_output = gsea_output, 
                exp.data.counts = pnoc008_counts_melt, 
                exp.data.tpm = pnoc008_tpm_melt,
                refData.counts = pbta_hgg_counts, 
                gene_set = dsigdb_geneset, 
                comparison = paste0("PBTA_HGG_", ncol(pbta_hgg_counts)))

# pnoc008 vs tcga gbm
gsea_output <- file.path(ref_dir, 'gsea', 'pnoc008_vs_tcga_gbm.rds')
gsea_enrichment(sample_to_add = sample_to_add, 
                gsea_output = gsea_output, 
                exp.data.counts = pnoc008_counts_melt, 
                exp.data.tpm = pnoc008_tpm_melt,
                refData.counts = tcga_gbm_counts, 
                gene_set = gene_set, 
                comparison = paste0("TCGA_GBM_", ncol(tcga_gbm_counts)))

# pbta comparisons only need to be run once
# pbta vs gtex brain
fname <- file.path(ref_dir, 'gsea', 'pbta_vs_gtex_brain.rds')
if(!file.exists(fname)){
  pbta_vs_gtex_brain <-  plyr::dlply(pbta_full_counts_melt, 
                                     .variables = "Sample", 
                                     .fun = function(x) run_rnaseq_analysis_edger(exp.data.counts = x, exp.data.tpm = pbta_full_tpm_melt, refData.counts = gtex_brain_counts, gene_set = gene_set, comparison = paste0("GTExBrain_", ncol(gtex_brain_counts))), .parallel = TRUE)
  saveRDS(pbta_vs_gtex_brain, file = fname)
}

# pbta vs pbta
fname <- file.path(ref_dir, 'gsea', 'pbta_vs_pbta.rds')
if(!file.exists(fname)){
  pbta_vs_pbta <- plyr::dlply(pbta_full_counts_melt, 
                              .variables = "Sample", 
                              .fun = function(x) run_rnaseq_analysis_edger(exp.data.counts = x, exp.data.tpm = pbta_full_tpm_melt, refData.counts = pbta_full_counts, gene_set = gene_set, comparison = paste0("PBTA_All_", ncol(pbta_full_counts))), .parallel = TRUE)
  saveRDS(pbta_vs_pbta, file = fname)
}

# pbta vs pbta hgg
fname <- file.path(ref_dir, 'gsea', 'pbta_vs_pbta_hgg.rds')
if(!file.exists(fname)){
  pbta_vs_pbta_hgg <- plyr::dlply(pbta_full_counts_melt, 
                                  .variables = "Sample", 
                                  .fun = function(x) run_rnaseq_analysis_edger(exp.data.counts = x, exp.data.tpm = pbta_full_tpm_melt, refData.counts = pbta_hgg_counts, gene_set = gene_set, comparison = paste0("PBTA_HGG_", ncol(pbta_hgg_counts))), .parallel = TRUE)
  saveRDS(pbta_vs_pbta_hgg, file = fname)
}

# tcga comparisons only need to be run once
# tcga gbm vs gtex brain
fname <- file.path(ref_dir, 'gsea', 'tcga_gbm_vs_gtex_brain.rds')
if(!file.exists(fname)){
  tcga_gbm_vs_gtex_brain <-  plyr::dlply(tcga_gbm_counts_melt, 
                                         .variables = "Sample", 
                                         .fun = function(x) run_rnaseq_analysis_edger(exp.data.counts = x, exp.data.tpm = tcga_gbm_tpm_melt, refData.counts = gtex_brain_counts, gene_set = gene_set, comparison = paste0("GTExBrain_", ncol(gtex_brain_counts))), .parallel = TRUE)
  saveRDS(tcga_gbm_vs_gtex_brain, file = fname)
}

# tcga gbm vs tcga gbm
fname <- file.path(ref_dir, 'gsea', 'tcga_gbm_vs_tcga_gbm.rds')
if(!file.exists(fname)){
  tcga_gbm_vs_tcga_gbm <-  plyr::dlply(tcga_gbm_counts_melt, 
                                       .variables = "Sample", 
                                       .fun = function(x) run_rnaseq_analysis_edger(exp.data.counts = x, exp.data.tpm = tcga_gbm_tpm_melt, refData.counts = tcga_gbm_counts, gene_set = gene_set, comparison = paste0("TCGA_GBM_", ncol(tcga_gbm_counts))), .parallel = TRUE)
  saveRDS(tcga_gbm_vs_tcga_gbm, file = fname)
}
