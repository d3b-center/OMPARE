# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'diffexpr_genes_barplot.R'))

# call function
fname <- file.path(topDir, "output", "diffexpr_genes_barplot_output.rds")
diffexpr_genes_barplot_output <- diffexpr_genes_barplot(rnaseq_analysis_output = rnaseq_analysis_output)
saveRDS(diffexpr_genes_barplot_output, file = fname)