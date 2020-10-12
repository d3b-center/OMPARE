# adult immune profile

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'immune_profile.R'))

# immune profile and save scores as well
adult_immune_profile <- immune_profile(fullmat = tcga.gbm.mat, 
                                           fname = file.path(topDir, 'output', 'immune_scores_adult.txt'))

# save output
saveRDS(adult_immune_profile, file = file.path(topDir, "output", "adult_immune_profile.rds"))