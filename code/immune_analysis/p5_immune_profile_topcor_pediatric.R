# pediatric immune profile (top correlated samples)

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "immune_analysis")
output_dir <- file.path(patient_dir, "output", "immune_analysis")
dir.create(output_dir, recursive = T, showWarnings = F)

# source functions
source(file.path(module_dir, "utils", 'immune_profile.R'))
source(file.path(module_dir, "utils", 'plot_immune_profile.R'))

# immune profile and save scores as well
fname <- file.path(output_dir, 'immune_scores_topcor_pediatric.rds')
if(!file.exists(fname)){
  pediatric_topcor_immune_profile <- immune_profile(fullmat = pbta_pnoc008_nn_tpm)
  saveRDS(pediatric_topcor_immune_profile, file = fname)
} else {
  pediatric_topcor_immune_profile <- readRDS(fname)
}

# plot immune scores
fname <- file.path(output_dir, 'immune_scores_topcor_pediatric.pdf')
if(!file.exists(fname)){
  pediatric_topcor_immune_profile <- plot_immune_profile(xcell_scores = pediatric_topcor_immune_profile)
  ggsave(filename = fname, plot = pediatric_topcor_immune_profile, 
         device = "pdf", width = 8, height = 10)
}