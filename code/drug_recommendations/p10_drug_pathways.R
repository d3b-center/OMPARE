suppressPackageStartupMessages({
  library(ggplot2)
  library(plyr)
  library(patchwork)
  library(dplyr)
})

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
module_dir <- file.path(root_dir, "code", "drug_recommendations")
output_dir <- file.path(patient_dir, "output", "drug_recommendations")
dir.create(output_dir, showWarnings = F, recursive = T)

# source functions
source(file.path(module_dir, "utils", "lincs_connectivity.R"))

# read edgeR output for all comparisons
input_dir <- file.path(patient_dir, 'output', 'rnaseq_analysis')

# patient vs normal tissue
patient_vs_normals_gsea <- readRDS(file.path(input_dir, 'patient_vs_normals_gsea.rds'))
patient_vs_normals_gsea$genes$comparison <- "patient_vs_normals"

# patient vs adult cancer
patient_vs_adult_gsea <- readRDS(file.path(input_dir, 'patient_vs_adult_gsea.rds'))
patient_vs_adult_gsea$genes$comparison <- "patient_vs_adult"

# patient vs pediatric cancer
patient_vs_pediatric_gsea <- readRDS(file.path(input_dir, 'patient_vs_pediatric_gsea.rds'))
patient_vs_pediatric_gsea$genes$comparison <- "patient_vs_pediatric"

# combine all outputs
genes_df <- rbind(patient_vs_normals_gsea$genes, patient_vs_adult_gsea$genes, patient_vs_pediatric_gsea$genes)

# LINCS-based similarity metric
drug_pathways_barplot <- plyr::dlply(.data = genes_df, 
                                     .variables = "comparison", 
                                     .fun = function(x) lincs_connectivity(input = x,
                                                                           method = "LINCS",
                                                                           trend_val = "down",
                                                                           output_dir = output_dir))
# save output
fname <- file.path(output_dir, "drug_pathways_barplot.pdf")
pdf(file = fname, width = 15, height = 10)
for(i in 1:length(drug_pathways_barplot)){
  print(drug_pathways_barplot[[i]])
}
dev.off()
