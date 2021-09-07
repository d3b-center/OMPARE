suppressPackageStartupMessages({
  library(ggplot2)
  library(plyr)
  library(patchwork)
})

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
module_dir <- file.path(root_dir, "code", "drug_recommendations")
output_dir <- file.path(patient_dir, "output", "drug_recommendations")
dir.create(output_dir, showWarnings = F, recursive = T)

# source functions
source(file.path(module_dir, "utils", "lincs_connectivity.R"))

# reference
gsea_dir <- file.path(data_dir, "gsea")

# read edgeR output for all comparisons
# gtex brain
pnoc008_vs_gtex_brain <- readRDS(file.path(gsea_dir, 'pnoc008_vs_gtex_brain.rds'))
pnoc008_vs_gtex_brain <- pnoc008_vs_gtex_brain[[patient]]

# pbta hgg
pnoc008_vs_pbta_hgg <- readRDS(file.path(gsea_dir, 'pnoc008_vs_pbta_hgg.rds'))
pnoc008_vs_pbta_hgg <- pnoc008_vs_pbta_hgg[[patient]]

# pbta all histologies
pnoc008_vs_pbta <- readRDS(file.path(gsea_dir, 'pnoc008_vs_pbta.rds'))
pnoc008_vs_pbta <- pnoc008_vs_pbta[[patient]]

# combine all outputs
genes_df <- rbind(pnoc008_vs_gtex_brain$genes, pnoc008_vs_pbta_hgg$genes, pnoc008_vs_pbta$genes)

# LINCS-based similarity metric
drug_pathways_barplot <- plyr::dlply(.data = genes_df, 
                                     .variables = "comparison", 
                                     .fun = function(x) lincs_connectivity(input = x,
                                                                           method = "LINCS",
                                                                           trend_val = "down",
                                                                           output_dir = output_dir))
fname <- file.path(output_dir, "drug_pathways_barplot.pdf")
ggsave(plot = wrap_plots(drug_pathways_barplot, ncol = 1), 
       filename = fname, 
       width = 23, height = 25, device = 'pdf')