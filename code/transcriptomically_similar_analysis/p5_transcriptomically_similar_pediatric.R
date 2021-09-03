# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "transcriptomically_similar_analysis")
output_dir <- file.path(patient_dir, "output", "transcriptomically_similar_analysis")
dir.create(output_dir, showWarnings = F, recursive = T)

# source functions
source(file.path(module_dir, "utils", 'transcriptomically_similar.R'))

# transcriptomically similar samples table
fname <- file.path(output_dir, "transciptomically_similar_pediatric.rds")
if(!file.exists(fname)){
  transciptomically_similar_pediatric <- transciptomically_similar(all_cor = pbta_pnoc008_nn_table, 
                                                                   clin_data = pbta_pnoc008_clinical)
  
  # save output
  saveRDS(transciptomically_similar_pediatric, file = fname)
} else {
  transciptomically_similar_pediatric <- readRDS(fname)
}