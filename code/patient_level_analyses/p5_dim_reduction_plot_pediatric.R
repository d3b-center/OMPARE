# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'dim_reduction_plot.R'))

# recurrent alterations
dim_reduction_plot_pediatric <- dim_reduction_plot(dat = pbta.embedding,
                   clindata = pbta.clinData,
                   study = "PBTA",
                   patient = sampleInfo$subjectID,
                   title =  "UMAP Correlation Clustering")
# save output
saveRDS(dim_reduction_plot_pediatric, file = file.path(topDir, "output", "dim_reduction_plot_pediatric.rds"))
