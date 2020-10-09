# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'tis_profile.R'))

# immune profile and save scores as well
tis_profile <- tis_profile(fname = file.path(topDir, 'output', 'tis_scores.txt'),  
                           score = "avg")

# save output
saveRDS(tis_profile, file = file.path(topDir, "output", "tis_profile.rds"))





