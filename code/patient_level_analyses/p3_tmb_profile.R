# TMB profile

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'tmb_profile.R')) 

# tmb profile
tmb_profile_output <- tmb_profile(pedTMBScores = pedTMB, 
                                  adultTMBScores = adultTMB, 
                                  TMB = tmb)

# save output
saveRDS(tmb_profile_output, file = file.path(topDir, "output", "tmb_profile_output.rds"))
