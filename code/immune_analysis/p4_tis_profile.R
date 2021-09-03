# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "immune_analysis")
output_dir <- file.path(patient_dir, "output", "immune_analysis")
dir.create(output_dir, recursive = T, showWarnings = F)

# source functions
source(file.path(module_dir, "utils", 'tis_profile.R'))
source(file.path(module_dir, "utils", 'plot_tis_profile.R'))

# save only scores
fname = file.path(output_dir, 'tis_scores.rds')
if(file.exists(fname)){
  tis_profile_output <- readRDS(fname)
} else {
  tis_profile_output <- tis_profile(patient_clinical = pnoc008_clinical, sampleInfo)
  saveRDS(tis_profile_output, file = fname)
}

# call function to plot data
fname <- file.path(output_dir, 'tis_scores.pdf')
if(!file.exists(fname)){
  tis_profile_output <- plot_tis_profile(tis_input = tis_profile_output, 
                                         score = "avg", 
                                         patient_of_interest = patient)
  ggsave(filename = fname, plot = tis_profile_output, 
         device = "pdf", width = 12, height = 8)
}

