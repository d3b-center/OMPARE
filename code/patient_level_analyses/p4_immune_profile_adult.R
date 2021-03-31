# adult immune profile

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'immune_profile.R'))
source(file.path(patient_level_analyses_utils, 'plot_immune_profile.R'))

# immune profile and save scores as well
fname <- file.path(topDir, 'output', 'immune_scores_adult.rds')
if(!file.exists(fname)){
  adult_immune_profile <- immune_profile(fullmat = tcga_gbm_pnoc008_immune_profile)
} else {
  adult_immune_profile <- readRDS(fname)
}

# plot immune scores
adult_immune_profile <- plot_immune_profile(xcell_scores = adult_immune_profile)