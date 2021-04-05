library(sva)
library(tidyverse)
library(ggridges)
library(ggplot2)
library(patchwork)
library(optparse)
library(dplyr)

# functions
source(file.path(patient_level_analyses_utils, 'quiet.R'))
source(file.path(patient_level_analyses_utils, 'batch_correct.R'))
source(file.path(patient_level_analyses_utils, "pubTheme.R"))

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))
pnoc008_dir <- file.path(ref_dir, 'pnoc008')
pbta_dir <- file.path(ref_dir, 'pbta')
gtex_dir <- file.path(ref_dir, 'gtex')

# output
dir <- file.path(topDir, 'output', 'drug_dge_density_plots')
dir.create(path = dir, showWarnings = F, recursive = T)

# create density plots and save output
dge_density_plots_helper <- function(combined_data, pnoc008_sample, dge_gene){
  vline <- combined_data %>%
    filter(sample == pnoc008_sample) %>%
    mutate(tpm = log2(tpm + 1)) %>%
    .$tpm
  combined_data <- combined_data %>%
    filter(study_id != "PNOC008")
  
  # compute stats
  d2 <- combined_data %>%
    group_by(study_id) %>%
    summarize(lower = quantile(log2(tpm + 1), probs = .025),
              median = quantile(log2(tpm + 1), probs = .5),
              upper = quantile(log2(tpm + 1), probs = .975))
  
  p <- ggplot(data = combined_data, aes(x = log2(tpm + 1))) +
    geom_density(aes(fill = study_id, alpha = 0.5)) +
    theme_Publication(base_size = 10) +
    facet_wrap(~study_id, scales = "free") +
    geom_vline(xintercept = vline, color = "red") + 
    geom_vline(data = d2, aes(xintercept = upper), color = "black", linetype = "dashed") + 
    geom_vline(data = d2, aes(xintercept = median), color = "black", linetype = "dashed") + 
    geom_vline(data = d2, aes(xintercept = lower), color = "black", linetype = "dashed") + 
    guides(alpha = F,  fill = F) +
    ggtitle(dge_gene) +
    xlab("log2 TPM") 
  
  return(p)
}

dge_density_plots <- function(topDir){
  # query file
  dge_all <- readRDS(file.path(topDir, 'output', 'transcriptome_drug_rec.rds'))
  dge_all <- dge_all %>% 
    dplyr::select(Gene, Comparison, logFC) %>%
    unique()
  dge_top <- dge_all %>%
    group_by(Comparison) %>%
    slice_max(logFC, n = 5) %>%
    mutate(top_dge = TRUE) %>%
    full_join(dge_all, by = c("Gene", "Comparison", "logFC"))
  dge_top$top_dge[is.na(dge_top$top_dge)] <- FALSE
  
  # pnoc008 clinical
  pnoc008_clinical <- readRDS(file.path(pnoc008_dir, 'pnoc008_clinical.rds'))
  pnoc008_clinical_remove <- pnoc008_clinical$Kids_First_Biospecimen_ID[pnoc008_clinical$subjectID != "PNOC008-29"]
  pnoc008_sample <- pnoc008_clinical$Kids_First_Biospecimen_ID[pnoc008_clinical$subjectID == "PNOC008-29"]
  tgen_samples <- c("BC", "BR", "BS", "FB", "XX")
  
  # reference file
  combined_tpm <- readRDS(file.path(pbta_dir, 'pbta-tgen-gtex-gene-expression-rsem-tpm-collapsed.combined.rds'))
  combined_tpm <- combined_tpm[,!colnames(combined_tpm) %in% c(tgen_samples, pnoc008_clinical_remove)]
  
  # histology file
  histology <- read_tsv(file.path(pbta_dir, 'pbta-histologies.tsv'))
  histology <- histology %>%
    filter(experimental_strategy == 'RNA-Seq',
           Kids_First_Biospecimen_ID %in% colnames(combined_tpm)) %>%
    mutate(study_id = ifelse(Kids_First_Biospecimen_ID %in% pnoc008_sample, "PNOC008", "PBTA")) %>%
    dplyr::select(Kids_First_Biospecimen_ID, RNA_library, short_histology, study_id) %>%
    unique()
  
  # gtex 
  gtex_clinical <- readRDS(file.path(gtex_dir, 'gtex_brain_clinical.rds'))
  gtex_clinical <- gtex_clinical %>%
    mutate(Kids_First_Biospecimen_ID = sample_id,
           RNA_library = library_name,
           short_histology = 'normal_brain',
           study = 'GTEx') %>%
    dplyr::select(Kids_First_Biospecimen_ID, RNA_library, short_histology, study_id)
  
  # combine both 
  combined_histology <- rbind(histology, gtex_clinical)
  combined_histology <- combined_histology %>%
    remove_rownames() %>%
    column_to_rownames('Kids_First_Biospecimen_ID')
  combined_histology <- combined_histology[colnames(combined_tpm),]
  
  # RNA_library and study_id constitute a batch
  combined_histology <- combined_histology %>%
    mutate(batch = paste0(RNA_library, "_", study_id))
  
  # batch correct 
  combined_tpm_corrected <- quiet(batch.correct(mat = combined_tpm, clin = combined_histology))
  
  # now subset to genes of interest
  combined_tpm_corrected <- combined_tpm_corrected[unique(dge_all$Gene),]
  combined_tpm_corrected <- melt(as.matrix(combined_tpm_corrected), varnames = c("gene", "sample"), value.name = "tpm")
  combined_tpm_corrected <- combined_tpm_corrected %>%
    inner_join(combined_histology %>%
                 rownames_to_column("sample"), by = "sample")
  
  # add pbta hgat 
  pbta_hgat <- combined_tpm_corrected %>%
    filter(short_histology  == "HGAT",
           study_id == "PBTA") %>%
    mutate(study_id = "PBTA_HGAT")
  combined_tpm_corrected <- rbind(combined_tpm_corrected, pbta_hgat)
  
  # full set of genes
  plist <- list()
  dge_genes = unique(dge_all$Gene)
  for(i in 1:length(dge_genes)){
    dge_gene = dge_genes[i]
    # for a combined plot with all studies
    mat <- combined_tpm_corrected %>%
      filter(gene == dge_gene)
    
    #  call function
    plist[[i]] <- dge_density_plots_helper(combined_data = mat, pnoc008_sample = pnoc008_sample, dge_gene = dge_gene)
    names(plist)[[i]] <- dge_gene
    
    # save output
    ggsave(plot = plist[[i]], filename = file.path(dir, paste0(dge_gene, '_drug_dge_density_plots.png')), width = 12, height = 5)
  }
  
  # output top genes as single png
  top_genes <- dge_top %>% filter(top_dge == TRUE) %>%
    .$Gene %>%
    unique
  top_genes <- plist[top_genes]
  ggsave(wrap_plots(top_genes, ncol = 3), 
         filename = file.path(dir, 'top_drug_dge_density_plots.png'), 
         width = 28, height = 14, device = 'png')
}

