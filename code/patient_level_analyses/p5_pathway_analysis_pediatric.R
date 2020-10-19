# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'pathway_analysis.R'))

# recurrent alterations
pathway_analysis_pediatric <- pathway_analysis(all_cor = pbta_pnoc008_nn_table, 
                                               prefix = "pediatric_", 
                                               comparison = "pediatric")

# save output
saveRDS(pathway_analysis_pediatric, file = file.path(topDir, "output", "pathway_analysis_pediatric.rds"))
