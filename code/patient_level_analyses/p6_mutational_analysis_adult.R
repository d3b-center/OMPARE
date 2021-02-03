# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'mutational_analysis.R'))

# recurrent alterations
mutational_analysis_adult <- mutational_analysis(top_cor = tcga_gbm_pnoc008_nn_tpm, 
                                                 key_clinical_findings_output = key_clinical_findings_output,
                                                 comparison = "adult")

# save output
saveRDS(mutational_analysis_adult, file = file.path(topDir, "output", "mutational_analysis_adult.rds"))
