# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'kaplan_meier.R'))

# recurrent alterations
kaplan_meier_adult <- kaplan_meier(all_cor = tcga.gbm.allCor, surv_data = tcga.gbm.survData)

# save output
saveRDS(kaplan_meier_adult, file = file.path(topDir, "output", "kaplan_meier_adult.rds"))
