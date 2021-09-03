# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "transcriptomically_similar_analysis")
output_dir <- file.path(patient_dir, "output", "transcriptomically_similar_analysis")
dir.create(output_dir, showWarnings = F, recursive = T)

# source functions
source(file.path(module_dir, "utils", 'pathway_analysis.R'))

# pathway enrichment top correlated samples
fname <- file.path(output_dir, "pathway_analysis_adult.rds")
if(!file.exists(fname)){
  pathway_analysis_adult <- pathway_analysis(all_cor = tcga_gbm_pnoc008_nn_table, 
                                             prefix = "adult_", 
                                             comparison = "adult",
                                             patient_of_interest = patient)
  
  # save output
  saveRDS(pathway_analysis_adult, file = fname)
} else {
  pathway_analysis_adult <- readRDS(fname)
}
