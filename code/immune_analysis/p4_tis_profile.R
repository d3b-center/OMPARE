# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "immune_analysis")
output_dir <- file.path(patient_dir, "output", "immune_analysis")
dir.create(output_dir, recursive = T, showWarnings = F)

# source functions
source(file.path(module_dir, "utils", 'tis_profile.R'))
source(file.path(module_dir, "utils", 'plot_tis_profile.R'))

# save only scores
fname <- file.path(output_dir, 'tis_scores.rds')
tis_profile_output <- tis_profile(pediatric_dir = "data/master_genomics", 
                                  adult_dir = "data/OpenPedCan-analysis/data",
                                  patient_of_interest = patient,
                                  norm_method = "quantile",
                                  collapse = FALSE)
saveRDS(tis_profile_output, file = fname)

# call function to plot data
fname <- file.path(output_dir, 'tis_scores.pdf')
tis_profile_output <- plot_tis_profile(tis_profile_output = tis_profile_output, 
                                       score = "avg")
ggsave(filename = fname, plot = tis_profile_output, device = "pdf", width = 12, height = 8)
