# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'mutational_analysis.R'))

# recurrent alterations
mutational_analysis_pediatric <- mutational_analysis(top_cor = pbta_pnoc008_nn_tpm, 
                                                     key_clinical_findings_output = key_clinical_findings_output,
                                                     clinical = pbta_pnoc008_clinical,
                                                     comparison = "pediatric")

# save output
saveRDS(mutational_analysis_pediatric, file = file.path(topDir, "output", "mutational_analysis_pediatric.rds"))
