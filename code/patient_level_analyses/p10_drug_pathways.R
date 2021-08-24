# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))
gsea_dir <- file.path(ref_dir, "gsea")

# source functions
source(file.path(patient_level_analyses_utils, 'lincs_connectivity.R'))

# patient of interest
pnoc008_patient <- gsub('.*/','',topDir)

# read edgeR output
# gtex brain
pnoc008_vs_gtex_brain <- readRDS(file.path(gsea_dir, 'pnoc008_vs_gtex_brain.rds'))
pnoc008_vs_gtex_brain <- pnoc008_vs_gtex_brain[[pnoc008_patient]]

# pbta hgg
pnoc008_vs_pbta_hgg <- readRDS(file.path(gsea_dir, 'pnoc008_vs_pbta_hgg.rds'))
pnoc008_vs_pbta_hgg <- pnoc008_vs_pbta_hgg[[pnoc008_patient]]

# pbta all histologies
pnoc008_vs_pbta <- readRDS(file.path(gsea_dir, 'pnoc008_vs_pbta.rds'))
pnoc008_vs_pbta <- pnoc008_vs_pbta[[pnoc008_patient]]

# combine all
genes_df <- rbind(pnoc008_vs_gtex_brain$genes, pnoc008_vs_pbta_hgg$genes, pnoc008_vs_pbta$genes)

# LINCS-based similarity metric
drug_pathways_barplot <- plyr::dlply(.data = genes_df, 
                                     .variables = "comparison", 
                                     .fun = function(x) lincs_connectivity(input = x,
                                                                           method = "LINCS",
                                                                           trend_val = "down",
                                                                           output_dir = file.path(topDir, "output")))
fname <- file.path(topDir, "output", "drug_pathways_barplot.pdf")
ggsave(plot = wrap_plots(drug_pathways_barplot, ncol = 1), 
       filename = fname, 
       width = 23, height = 25, device = 'pdf')
