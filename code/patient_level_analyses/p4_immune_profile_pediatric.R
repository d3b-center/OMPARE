# pediatric immune profile

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'immune_profile.R'))

# immune profile and save scores as well
pediatric_immune_profile <- immune_profile(fullmat = pbta.mat.all, 
                                           fname = file.path(topDir, 'output', 'immune_scores_pediatric.txt'))

# save output
saveRDS(pediatric_immune_profile, file = file.path(topDir, "output", "pediatric_immune_profile.rds"))
