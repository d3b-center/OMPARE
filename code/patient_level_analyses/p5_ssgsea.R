# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'ssgsea.R'))
source(file.path(patient_level_analyses_utils, 'plot_ssgsea.R'))

# recurrent alterations
fname <- file.path(topDir, 'output', 'ssgsea_scores_pediatric.rds')
if(!file.exists(fname)){
  ssgsea_pediatric <- ssgsea(top_cor = pbta_pnoc008_nn_tpm)
  saveRDS(ssgsea_pediatric, file = fname)
} else {
  ssgsea_pediatric <- readRDS(fname)
}

# plot ssgea
fname <- file.path(topDir, 'output', 'ssgsea_scores_pediatric.pdf')
if(!file.exists(fname)){
  ssgsea_pediatric <- plot_ssgsea(ssgsea_pediatric)
  ggsave(filename = fname, plot = ssgsea_pediatric, 
         device = "pdf", width = 15, height = 14)
}
