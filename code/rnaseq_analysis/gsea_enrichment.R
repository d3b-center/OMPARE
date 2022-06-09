# Author: Komal S. Rathi
# function to do gsea between patient of interest vs normal tissue, adult cancer and pediatric cancer

suppressPackageStartupMessages({
  library(tidyverse)
  library(GSEABase)
  library(reshape2)
})

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
module_dir <- file.path(root_dir, "code", "rnaseq_analysis")

gsea_enrichment <- function(normal_tissue, adult_cancer, pediatric_cancer, output_dir, tpm_data, count_data){
  
  # format input tpm data
  tpm_data_melted <- tpm_data %>% 
    rownames_to_column("gene_symbol") %>%
    gather(sample, tpm, -c('gene_symbol'))
  
  # format input count data 
  count_data_melted <- count_data %>% 
    rownames_to_column("gene_symbol") %>%
    gather(sample, expected_count, -c('gene_symbol'))
  
  # patient of interest
  patient_of_interest <- unique(tpm_data_melted$sample)
  
  # source function for RNA-seq diffexpr & pathway analysis
  source(file.path(module_dir, "utils", "rnaseq_analysis_edgeR.R"))
  
  # cancer Genes
  cancer_genes <- readRDS(file.path(data_dir, 'cancer_gene_list.rds'))
  
  # genesets 
  gene_set <- getGmt(file.path(data_dir, 'msigdb', 'c2.cp.reactome.v6.0.symbols.gmt'), collectionType = BroadCollection(), geneIdType = SymbolIdentifier())
  gene_set <- geneIds(gene_set)
  
  # Dataset1: GTex normals
  # counts
  normal_tissue_dir <- file.path(root_dir, "data", "normal_data")
  normal_counts <- list.files(path = normal_tissue_dir, pattern = "*counts.rds", full.names = T)
  normal_counts <- readRDS(normal_counts)
  
  # Dataset2: Adult cancer
  # counts
  adult_cancer_dir <- file.path(root_dir, "data", "adult_data")
  adult_counts <- list.files(path = adult_cancer_dir, pattern = "*counts.rds", full.names = T)
  adult_counts <- readRDS(adult_counts)
  
  # Dataset3: Pediatric cancer
  # counts
  pediatric_cancer_dir <- file.path(root_dir, "data", "pediatric_data")
  pediatric_counts <- list.files(path = pediatric_cancer_dir, pattern = "*counts.rds", full.names = T)
  pediatric_counts <- readRDS(pediatric_counts)
  
  # remove patient of interest
  pediatric_counts <- pediatric_counts %>%
    dplyr::select(-c(patient_of_interest))
  
  # patient vs normal data
  gsea_vs_normal <- file.path(output_dir, 'patient_vs_normals_gsea.rds')
  if(!file.exists(gsea_vs_normal)){
    res_gsea_vs_normal <- run_rnaseq_analysis_edger(exp.data.counts = count_data_melted, 
                                                    exp.data.tpm = tpm_data_melted, 
                                                    refData.counts = normal_counts, 
                                                    gene_set = gene_set, 
                                                    comparison = paste0("Normal: ", normal_tissue),
                                                    cancer_genes = cancer_genes,
                                                    single_sample = FALSE,
                                                    sample_name = patient_of_interest)
    saveRDS(res_gsea_vs_normal, file = gsea_vs_normal)
    
    # function to merge degs for patient vs normal data
    deg_vs_normal <- res_gsea_vs_normal$genes
    deg_vs_normal <- deg_vs_normal %>%
      dplyr::select(genes, diff_expr) %>%
      mutate(sample_name = patient_of_interest) 
    saveRDS(deg_vs_normal, file = file.path(output_dir, "patient_vs_normal_degs.rds"))
  } else {
    res_gsea_vs_normal <- readRDS(gsea_vs_normal)
  }
  
  # patient vs pediatric
  gsea_vs_pediatric <- file.path(output_dir, 'patient_vs_pediatric_gsea.rds')
  if(!file.exists(gsea_vs_pediatric)){
    res_gsea_vs_pediatric <- run_rnaseq_analysis_edger(exp.data.counts = count_data_melted, 
                                                       exp.data.tpm = tpm_data_melted, 
                                                       refData.counts = pediatric_counts, 
                                                       gene_set = gene_set, 
                                                       comparison = paste0("Pediatric: ", pediatric_cancer),
                                                       cancer_genes = cancer_genes,
                                                       single_sample = FALSE,
                                                       sample_name = patient_of_interest)
    saveRDS(res_gsea_vs_pediatric, file = gsea_vs_pediatric)
  } else {
    res_gsea_vs_pediatric <- readRDS(gsea_vs_pediatric)
  }
  
  # patient vs adult
  gsea_vs_adult <- file.path(output_dir, 'patient_vs_adult_gsea.rds')
  if(!file.exists(gsea_vs_adult)){
    res_gsea_vs_adult <- run_rnaseq_analysis_edger(exp.data.counts = count_data_melted, 
                                                   exp.data.tpm = tpm_data_melted, 
                                                   refData.counts = adult_counts, 
                                                   gene_set = gene_set, 
                                                   comparison = paste0("Adult: ", adult_cancer),
                                                   cancer_genes = cancer_genes,
                                                   single_sample = FALSE,
                                                   sample_name = patient_of_interest)
    saveRDS(res_gsea_vs_adult, file = gsea_vs_adult)
  } else {
    res_gsea_vs_adult <- readRDS(gsea_vs_adult)
  }
  
  # write out to text files
  # up/down pathways
  pathway_df <- rbind(res_gsea_vs_normal$pathways, 
                      res_gsea_vs_pediatric$pathways, 
                      res_gsea_vs_adult$pathways)
  pathway_df <- pathway_df %>%
    group_by(pathway, direction) %>%
    mutate(Freq = n()) %>%
    as.data.frame()
  pathway_df_up <- pathway_df %>%
    filter(direction == "up") %>%
    as.data.frame() %>%
    write_tsv(file = file.path(output_dir, "pathways_up.txt"))
  pathway_df_down <- pathway_df %>%
    filter(direction == "down") %>%
    as.data.frame() %>%
    write_tsv(file = file.path(output_dir, "pathways_down.txt"))
  
  # up/down genes
  genes_df <- rbind(res_gsea_vs_normal$genes, 
                    res_gsea_vs_pediatric$genes, 
                    res_gsea_vs_adult$genes)
  genes_df <- genes_df %>%
    group_by(genes, diff_expr) %>%
    mutate(Freq = n()) %>%
    as.data.frame()
  genes_df_up <- genes_df %>%
    filter(diff_expr == "up") %>%
    as.data.frame() %>%
    write_tsv(file = file.path(output_dir, "genes_up.txt"))
  genes_df_down <- genes_df %>%
    filter(diff_expr == "down") %>%
    as.data.frame() %>%
    write_tsv(file = file.path(output_dir, "genes_down.txt"))
}

