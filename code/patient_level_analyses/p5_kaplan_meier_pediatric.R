# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'kaplan_meier.R'))

# recurrent alterations
kaplan_meier_pediatric <- kaplan_meier(all_cor = pbta_hgat_pnoc008_nn_table, 
                                       surv_data = pbta_survival)

# save output
saveRDS(kaplan_meier_pediatric, file = file.path(topDir, "output", "kaplan_meier_pediatric.rds"))
