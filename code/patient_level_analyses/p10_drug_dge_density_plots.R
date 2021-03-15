# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'drug_dge_density_plots.R'))

# call function
drug_rec_dge_density_plots <- dge_density_plots(topDir = topDir)
