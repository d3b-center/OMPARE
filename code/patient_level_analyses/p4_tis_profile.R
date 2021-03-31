# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'tis_profile.R'))
source(file.path(patient_level_analyses_utils, 'plot_tis_profile.R'))

# save only scores
fname = file.path(topDir, 'output', 'tis_scores.rds')
if(file.exists(fname)){
  tis_profile_output <- readRDS(fname)
} else {
  tis_profile_output <- tis_profile(patient_clinical = pnoc008_clinical, sampleInfo)
  saveRDS(tis_profile_output, file = fname)
}

# call function to plot data
tis_profile_output <- plot_tis_profile(tis_input = tis_profile_output, score = "avg", sampleInfo = sampleInfo)
