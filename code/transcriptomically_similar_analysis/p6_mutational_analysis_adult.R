# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "transcriptomically_similar_analysis")
output_dir <- file.path(patient_dir, "output", "transcriptomically_similar_analysis")
dir.create(output_dir, showWarnings = F, recursive = T)

# source functions
source(file.path(module_dir, "utils", 'mutational_analysis.R'))

# recurrent alterations
fname <- file.path(output_dir, "mutational_analysis_adult.rds")
if(!file.exists(fname)){
  mutational_analysis_adult <- mutational_analysis(top_cor = tcga_gbm_pnoc008_nn_tpm, 
                                                   key_clinical_findings_output = key_clinical_findings_output,
                                                   clinical = tcga_gbm_pnoc008_clinical,
                                                   comparison = "adult")
  
  # save output
  saveRDS(mutational_analysis_adult, file = fname)
} else {
  mutational_analysis_adult <- readRDS(fname)
}