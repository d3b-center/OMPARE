# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "transcriptomically_similar_analysis")
output_dir <- file.path(patient_dir, "output", "transcriptomically_similar_analysis")
dir.create(output_dir, showWarnings = F, recursive = T)

# source functions
source(file.path(module_dir, "utils", "dim_reduction_plot.R"))

# umap clustering plot
fname <- file.path(output_dir, "dim_reduction_plot_adult.rds")
dim_reduction_plot_adult <- dim_reduction_plot(dat = tcga_gbm_pnoc008_embedding,
                                               clindata = tcga_gbm_pnoc008_clinical,
                                               study = "TCGA",
                                               patient = patient,
                                               title =  "UMAP Correlation Clustering")
# save output
saveRDS(dim_reduction_plot_adult, file = fname)