# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'dim_reduction_plot.R'))

# recurrent alterations
dim_reduction_plot_adult <- dim_reduction_plot(dat = tcga_gbm_pnoc008_embedding,
                                               clindata = tcga_gbm_pnoc008_clinical,
                                               study = "TCGA",
                                               patient = sampleInfo$subjectID,
                                               title =  "UMAP Correlation Clustering")
# save output
saveRDS(dim_reduction_plot_adult, file = file.path(topDir, "output", "dim_reduction_plot_adult.rds"))
