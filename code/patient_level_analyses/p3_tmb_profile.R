# TMB profile

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'tmb_profile.R')) 

# tmb profile
tmb_profile_output <- tmb_profile(pbta_tmb = pbta_tmb, 
                                  tcga_in_pbta_tmb = tcga_in_pbta_tmb, 
                                  tcga_not_in_pbta_tmb = tcga_not_in_pbta_tmb,
                                  tmb_bed_file = tmb_bed_file,
                                  TMB = tmb)

# save output
saveRDS(tmb_profile_output, file = file.path(topDir, "output", "tmb_profile_output.rds"))
