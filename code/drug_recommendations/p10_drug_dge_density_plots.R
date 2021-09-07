suppressPackageStartupMessages({
  library(sva)
  library(reshape2)
  library(tidyverse)
  library(ggplot2)
  library(dplyr)
})

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
pnoc008_dir <- file.path(data_dir, "pnoc008")
pbta_dir <- file.path(data_dir, "pbta")
gtex_dir <- file.path(data_dir, "gtex")

# output directory
module_dir <- file.path(root_dir, "code", "drug_recommendations")
output_dir <- file.path(patient_dir, "output", "drug_recommendations")
output_dir <- file.path(output_dir, "drug_dge_density_plots")
dir.create(output_dir, showWarnings = F, recursive = T)

# source functions
source(file.path(module_dir, "utils", "drug_dge_density_plots.R"))
source(file.path(root_dir, "code", "utils", "quiet.R"))
source(file.path(root_dir, "code", "utils", "batch_correct.R"))

# call function
fname <- file.path(output_dir, "top_drug_dge_density_plots.pdf")
if(!file.exists(fname)){
  
  # query file
  transcriptomic_drug_rec <- readRDS(file.path(patient_dir, "output", "drug_recommendations", "transcriptome_drug_rec.rds"))
  dge_all <- transcriptomic_drug_rec %>% 
    dplyr::select(Gene, Comparison, logFC) %>%
    unique()
  dge_top <- dge_all %>%
    group_by(Comparison) %>%
    slice_max(logFC, n = 5) %>%
    dplyr::mutate(top_dge = TRUE) %>%
    full_join(dge_all, by = c("Gene", "Comparison", "logFC"))
  dge_top$top_dge[is.na(dge_top$top_dge)] <- FALSE
  
  # pnoc008 clinical
  pnoc008_clinical <- readRDS(file.path(root_dir, "data", "pnoc008", "pnoc008_clinical.rds"))
  subject_id <- pnoc008_clinical %>% 
    filter(subjectID == patient) %>%
    .$subjectID
  sample_of_interest <- pnoc008_clinical %>% 
    filter(subjectID == patient) %>%
    mutate(RNA_library = library_name,
           short_histology = "HGAT") %>%
    dplyr::select(Kids_First_Biospecimen_ID, RNA_library, short_histology, study_id)
  
  # pbta
  # pbta expression matrix (n = 1035)
  pbta_tpm <- readRDS(file.path(pbta_dir, "pbta-gene-expression-rsem-tpm-collapsed.polya.stranded.rds"))
  
  # pbta histology file
  histology <- read_tsv(file.path(pbta_dir, 'pbta-histologies.tsv'))
  histology <- histology %>%
    filter(experimental_strategy == 'RNA-Seq',
           Kids_First_Biospecimen_ID %in% colnames(pbta_tpm)) %>%
    mutate(study_id = "PBTA") %>%
    dplyr::select(Kids_First_Biospecimen_ID, RNA_library, short_histology, study_id) %>%
    unique()
  
  # gtex brain tpm matrix
  gtex_tpm <- readRDS(file.path(gtex_dir, "gtex_brain_tpm.rds"))
  
  # gtex clinical file
  gtex_clinical <- readRDS(file.path(gtex_dir, "gtex_brain_clinical.rds"))
  gtex_clinical <- gtex_clinical %>%
    mutate(Kids_First_Biospecimen_ID = sample_id,
           RNA_library = library_name,
           short_histology = 'normal_brain',
           study = 'GTEx') %>%
    dplyr::select(Kids_First_Biospecimen_ID, RNA_library, short_histology, study_id)
  
  # combine clinical files 
  combined_histology <- rbind(histology, gtex_clinical, sample_of_interest)
  combined_histology <- combined_histology %>%
    remove_rownames() %>%
    column_to_rownames('Kids_First_Biospecimen_ID')
  
  # combine tpm matrices
  tpm_data <- expData %>%
    column_to_rownames('gene_symbol')
  colnames(tpm_data)[1] <- sample_of_interest$Kids_First_Biospecimen_ID
  common_genes = intersect(rownames(pbta_tpm), intersect(rownames(gtex_tpm), rownames(tpm_data)))
  combined_tpm <- cbind(pbta_tpm[common_genes,], gtex_tpm[common_genes,], tpm_data[common_genes,sample_of_interest$Kids_First_Biospecimen_ID, drop = F])
  combined_tpm <- combined_tpm[,colnames(combined_tpm) %in% rownames(combined_histology)]
  combined_histology <- combined_histology[colnames(combined_tpm),]
  
  # RNA_library and study_id constitute a batch
  combined_histology <- combined_histology %>%
    mutate(batch = paste0(RNA_library, "_", study_id))
  
  # batch correct 
  combined_tpm_corrected <- quiet(batch_correct(mat = combined_tpm, clin = combined_histology))
  
  # now subset to genes of interest
  combined_tpm_corrected <- combined_tpm_corrected[unique(dge_all$Gene),]
  combined_tpm_corrected <- melt(as.matrix(combined_tpm_corrected), varnames = c("gene", "sample"), value.name = "tpm")
  combined_tpm_corrected <- combined_tpm_corrected %>%
    inner_join(combined_histology %>%
                 rownames_to_column("sample"), by = "sample")
  
  # add pbta hgat (n = 189)
  pbta_hgat <- combined_tpm_corrected %>%
    filter(short_histology  == "HGAT",
           study_id == "PBTA") %>%
    dplyr::mutate(study_id = "PBTA_HGAT")
  combined_tpm_corrected <- rbind(combined_tpm_corrected, pbta_hgat)
  
  # create plots for full set of genes
  plist <- list()
  dge_genes = unique(dge_all$Gene)
  for(i in 1:length(dge_genes)){
    dge_gene = dge_genes[i]
    # for a combined plot with all studies
    mat <- combined_tpm_corrected %>%
      filter(gene == dge_gene)
    
    #  call function
    plist[[i]] <- drug_dge_density_plots(combined_data = mat, pnoc008_sample = sample_of_interest$Kids_First_Biospecimen_ID, dge_gene = dge_gene)
    names(plist)[[i]] <- dge_gene
    
    # save output
    ggsave(plot = plist[[i]], filename = file.path(output_dir, paste0(dge_gene, '_drug_dge_density_plots.png')), width = 12, height = 5)
  }
  
  # output top genes as single pdf
  top_genes <- dge_top %>% filter(top_dge == TRUE) %>%
    .$Gene %>%
    unique
  pdf(file = fname, width = 12, height = 5)
  for(i in 1:length(top_genes)){
    print(plist[top_genes[i]])
  }
  dev.off()
}
