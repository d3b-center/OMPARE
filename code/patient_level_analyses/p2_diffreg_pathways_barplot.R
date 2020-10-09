# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'diffreg_pathways_barplot.R'))

# call function
fname <- file.path(topDir, "output", "diffreg_pathways_barplot_output.rds")
diffreg_pathways_barplot_output <- diffreg_pathways_barplot(rnaseq_analysis_output = rnaseq_analysis_output)
saveRDS(diffreg_pathways_barplot_output, file = fname)
