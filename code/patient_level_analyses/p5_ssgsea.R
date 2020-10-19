# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'ssgsea.R'))

# recurrent alterations
ssgsea_pediatric <- ssgsea(top_cor = pbta_pnoc008_nn_tpm, 
                           fname = file.path(topDir, 'output', 'ssgsea_scores_pediatric.txt'))

# save output
saveRDS(ssgsea_pediatric, file = file.path(topDir, "output", "ssgsea_pediatric.rds"))
