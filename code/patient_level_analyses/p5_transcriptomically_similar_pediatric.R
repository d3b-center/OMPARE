# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'transcriptomically_similar.R'))

# recurrent alterations
transciptomically_similar_pediatric <- transciptomically_similar(all_cor = pbta_pnoc008_nn_table, 
                                                                 clin_data = pbta_pnoc008_clinical)

# save output
saveRDS(transciptomically_similar_pediatric, file = file.path(topDir, "output", "transciptomically_similar_pediatric.rds"))
