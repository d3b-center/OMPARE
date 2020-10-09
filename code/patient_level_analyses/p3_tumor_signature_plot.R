# TMB profile

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'tumor_signature_plot.R')) 

# tmb profile
tumor_signature_output <- tumor_signature_plot()

# save output
saveRDS(tumor_signature_output, file = file.path(topDir, "output", "tumor_signature_output.rds"))
