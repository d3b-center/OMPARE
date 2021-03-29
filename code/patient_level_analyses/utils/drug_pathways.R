# drug pathways
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(dplyr))

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))
dsigdb_dir <- file.path(ref_dir, "dsigdb")

# source functions
source(file.path(patient_level_analyses_utils, 'diffreg_pathways_barplot.R'))

drug_pathways <- function(pnoc008_patient, output_dir){
  # read output from gsea enrichment
  # gtex brain
  pnoc008_vs_gtex_brain <- readRDS(file.path(dsigdb_dir, 'pnoc008_vs_gtex_brain.rds'))
  pnoc008_vs_gtex_brain <- pnoc008_vs_gtex_brain[[pnoc008_patient]]
  
  # pbta hgg
  pnoc008_vs_pbta_hgg <- readRDS(file.path(dsigdb_dir, 'pnoc008_vs_pbta_hgg.rds'))
  pnoc008_vs_pbta_hgg <- pnoc008_vs_pbta_hgg[[pnoc008_patient]]
  
  # pbta all histologies
  pnoc008_vs_pbta <- readRDS(file.path(dsigdb_dir, 'pnoc008_vs_pbta.rds'))
  pnoc008_vs_pbta <- pnoc008_vs_pbta[[pnoc008_patient]]
  
  # up/down pathways
  pathway_df <- rbind(pnoc008_vs_gtex_brain$pathways, pnoc008_vs_pbta_hgg$pathways, pnoc008_vs_pbta$pathways)
  pathway_df <- pathway_df %>%
    mutate(genes = as.character(genes)) %>%
    group_by(pathway, direction) %>%
    mutate(Freq = n()) %>%
    as.data.frame()
  pathway_df.up <- pathway_df %>%
    filter(direction == "up") %>%
    as.data.frame()
  pathway_df.down <- pathway_df %>%
    filter(direction == "down") %>%
    as.data.frame()
  
  # up/down genes
  genes_df <- rbind(pnoc008_vs_gtex_brain$genes, pnoc008_vs_pbta_hgg$genes, pnoc008_vs_pbta$genes)
  genes_df <- genes_df %>%
    group_by(genes, diff_expr) %>%
    mutate(Freq = n()) %>%
    as.data.frame()
  genes_df.up <- genes_df %>%
    filter(diff_expr == "up") %>%
    as.data.frame()
  genes_df.down <- genes_df %>%
    filter(diff_expr == "down") %>%
    as.data.frame()
  
  # write output
  write.table(x = pathway_df.up, file = file.path(output_dir, "dsigdb_pathways_up.txt"), row.names = F, quote = F, sep = "\t")
  write.table(x = genes_df.up, file = file.path(output_dir, "dsigdb_de_genes_up.txt"), row.names = F, quote = F, sep = "\t")
  write.table(x = pathway_df.down, file = file.path(output_dir, "dsigdb_pathways_down.txt"), row.names = F, quote = F, sep = "\t")
  write.table(x = genes_df.down, file = file.path(output_dir, "dsigdb_de_genes_down.txt"), row.names = F, quote = F, sep = "\t")
  
  # barplot
  diffreg_pathways_barplot_gtex <- diffreg_pathways_barplot(pathway_df.up, pathway_df.down, comparison_study = 'GTExBrain_1152')
  diffreg_pathways_barplot_pbta_hgg <- diffreg_pathways_barplot(pathway_df.up, pathway_df.down, comparison_study = 'PBTA_HGG_182')
  diffreg_pathways_barplot_pbta <- diffreg_pathways_barplot(pathway_df.up, pathway_df.down, comparison_study = 'PBTA_All_1028')
  diffreg_pathways_barplot_output <- list(diffreg_pathways_barplot_gtex = diffreg_pathways_barplot_gtex,
                                          diffreg_pathways_barplot_pbta_hgg = diffreg_pathways_barplot_pbta_hgg,
                                          diffreg_pathways_barplot_pbta = diffreg_pathways_barplot_pbta)
  
  return(diffreg_pathways_barplot_output)
}
