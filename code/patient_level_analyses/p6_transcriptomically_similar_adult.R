# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'transcriptomically_similar.R'))

# recurrent alterations
transciptomically_similar_adult <- transciptomically_similar(all_cor = tcga_gbm_pnoc008_nn_table, 
                                                             clin_data = tcga_gbm_pnoc008_clinical)

# save output
saveRDS(transciptomically_similar_adult, file = file.path(topDir, "output", "transciptomically_similar_adult.rds"))
