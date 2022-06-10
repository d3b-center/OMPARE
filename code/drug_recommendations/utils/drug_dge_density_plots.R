# function to combine TPM matrices and generate drug dge density plots 
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggridges)
  library(ggplot2)
  library(dplyr)
})

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "quiet.R"))
source(file.path(root_dir, "code", "utils", "batch_correct.R"))
source(file.path(root_dir, "code", "utils", "pubTheme.R"))

drug_dge_density_plots <- function(patient_of_interest, dge_all, output_dir){
  
  # tpm for pediatric comparator
  pediatric_cancer_dir <- file.path(root_dir, "data", "pediatric_data")
  pediatric_tpm <- list.files(path = pediatric_cancer_dir, pattern = "tpm.rds", full.names = T)
  pediatric_tpm <- readRDS(pediatric_tpm)
  pediatric_clin <- list.files(path = pediatric_cancer_dir, pattern = "histologies.tsv", full.names = T)
  pediatric_clin <- read.delim(pediatric_clin)
  pediatric_clin <- pediatric_clin %>%
    filter(Kids_First_Biospecimen_ID %in% colnames(pediatric_tpm)) %>%
    mutate(type = "Pediatric tumors")
  
  # tpm for adult comparator
  adult_cancer_dir <- file.path(root_dir, "data", "adult_data")
  adult_tpm <- list.files(path = adult_cancer_dir, pattern = "tpm.rds", full.names = T)
  adult_tpm <- readRDS(adult_tpm)
  adult_clin <- list.files(path = adult_cancer_dir, pattern = "histologies.tsv", full.names = T)
  adult_clin <- read.delim(adult_clin)
  adult_clin <- adult_clin %>%
    filter(Kids_First_Biospecimen_ID %in% colnames(adult_tpm)) %>%
    mutate(type = "Adult tumors")
  
  # tpm for normal tissues
  normal_tissue_dir <- file.path(root_dir, "data", "normal_data")
  normal_tpm <- list.files(path = normal_tissue_dir, pattern = "tpm.rds", full.names = T)
  normal_tpm <- readRDS(normal_tpm)
  normal_clin <- list.files(path = normal_tissue_dir, pattern = "histologies.tsv", full.names = T)
  normal_clin <- read.delim(normal_clin)
  normal_clin <- normal_clin %>%
    filter(Kids_First_Biospecimen_ID %in% colnames(normal_tpm)) %>%
    mutate(type = "Normal tissues")
  
  
  # combine tpm for all data
  combined_tpm <- pediatric_tpm %>%
    rownames_to_column("gene_symbol") %>%
    inner_join(adult_tpm %>%
                 rownames_to_column("gene_symbol"), by = "gene_symbol") %>%
    inner_join(normal_tpm %>%
                 rownames_to_column("gene_symbol"), by = "gene_symbol") %>%
    column_to_rownames("gene_symbol")
  
  # combine histologies for all data
  common_cols <- intersect(intersect(colnames(pediatric_clin), colnames(adult_clin)), colnames(normal_clin))
  combined_histology <- rbind(pediatric_clin %>% select(common_cols), 
                         adult_clin %>% select(common_cols), 
                         normal_clin %>% select(common_cols)) %>%
    column_to_rownames('Kids_First_Biospecimen_ID') %>%
    mutate(batch = paste0(RNA_library, "_", cohort))
  
  # batch correct 
  combined_tpm <- combined_tpm[,rownames(combined_histology)]
  combined_tpm_corrected <- quiet(batch_correct(mat = combined_tpm, clin = combined_histology))
  
  # subset the corrected matrix to genes of interest
  combined_tpm_corrected <- combined_tpm_corrected[unique(dge_all$Gene),]
  combined_tpm_corrected <- melt(as.matrix(combined_tpm_corrected), varnames = c("gene", "sample"), value.name = "tpm")
  combined_tpm_corrected <- combined_tpm_corrected %>%
    inner_join(combined_histology %>%
                 rownames_to_column("sample"), by = "sample")
  
  # create plots for full set of genes
  plist <- list()
  dge_genes <- unique(dge_all$Gene)
  for(i in 1:length(dge_genes)){
    
    # gene of interest
    dge_gene <- dge_genes[i]
    
    # for a combined plot with all studies
    mat <- combined_tpm_corrected %>%
      filter(gene == dge_gene)
    
    # patient of interest
    vline <- mat %>%
      filter(sample == patient_of_interest) %>%
      dplyr::mutate(tpm = log2(tpm + 1)) %>%
      .$tpm
    
    # matrix of all data except patient of interest
    mat <- mat %>%
      filter(!sample %in% patient_of_interest)
    
    # add labels
    mat <- mat %>%
      group_by(type, short_histology) %>%
      dplyr::mutate(group = paste0(type, ": ", short_histology, "\n(n = ",n(),")"))
    
    # compute stats
    stats <- mat %>%
      group_by(group) %>%
      summarize(lower = quantile(log2(tpm + 1), probs = .025),
                median = quantile(log2(tpm + 1), probs = .5),
                upper = quantile(log2(tpm + 1), probs = .975))
    
    # create density plot and add to a list
    plist[[i]] <- ggplot(data = mat, aes(x = log2(tpm + 1))) +
      geom_density(aes(fill = group, alpha = 0.5)) +
      theme_Publication(base_size = 10) +
      facet_wrap(~group, scales = "free") +
      geom_vline(xintercept = vline, color = "red") + 
      geom_vline(data = stats, aes(xintercept = upper), color = "black", linetype = "dashed") + 
      geom_vline(data = stats, aes(xintercept = median), color = "black", linetype = "dashed") + 
      geom_vline(data = stats, aes(xintercept = lower), color = "black", linetype = "dashed") + 
      guides(alpha = "none",  fill = "none") +
      ggtitle(dge_gene) +
      xlab("log2 TPM") 
    names(plist)[[i]] <- dge_gene
    
    # save output
    ggsave(plot = plist[[i]], filename = file.path(output_dir, paste0(dge_gene, '_drug_dge_density_plots.png')), width = 12, height = 5)
  }
  
  # filter to top genes
  top_genes <- dge_all %>%
    group_by(Comparison) %>%
    slice_max(logFC, n = 5) %>%
    dplyr::mutate(top_dge = TRUE) %>%
    filter(!is.na(top_dge)) %>%
    .$Gene %>%
    unique()
  
  fname <- file.path(output_dir, "top_drug_dge_density_plots.pdf")
  pdf(file = fname, width = 12, height = 5)
  for(i in 1:length(top_genes)){
    print(plist[top_genes[i]])
  }
  dev.off()
}

