# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'pathway_analysis.R'))

# pathway enrichment top correlated samples
pathway_analysis_adult <- pathway_analysis(all_cor = tcga_gbm_pnoc008_nn_table, 
                                           prefix = "adult_", comparison = "adult")

# save output
saveRDS(pathway_analysis_adult, file = file.path(topDir, "output", "pathway_analysis_adult.rds"))
