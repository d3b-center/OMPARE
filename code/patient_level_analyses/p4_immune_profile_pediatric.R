# pediatric immune profile

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'immune_profile.R'))
source(file.path(patient_level_analyses_utils, 'plot_immune_profile.R'))

# immune profile and save scores as well
fname <- file.path(topDir, 'output', 'immune_scores_pediatric.rds')
if(!file.exists(fname)){
  pediatric_immune_profile <- immune_profile(fullmat = pbta_pnoc008_immune_profile)
  saveRDS(pediatric_immune_profile, file = fname)
} else {
  pediatric_immune_profile <- readRDS(fname)
}

# plot immune scores
fname <- file.path(topDir, 'output', 'immune_scores_pediatric.pdf')
if(!file.exists(fname)){
  pediatric_immune_profile <- plot_immune_profile(xcell_scores = pediatric_immune_profile)
  ggsave(filename = fname, plot = pediatric_immune_profile, 
         device = "pdf", width = 8, height = 10)
}
