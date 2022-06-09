# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "transcriptomically_similar_analysis")
output_dir <- file.path(patient_dir, "output", "transcriptomically_similar_analysis")
dir.create(output_dir, showWarnings = F, recursive = T)

# source functions
source(file.path(module_dir, "utils", 'transcriptomically_similar.R'))

# load inputs
nn_table <- file.path(output_dir, "pediatric_all_nn_table.rds")
nn_table <- readRDS(nn_table)
pediatric_patient_clinical <-  file.path(output_dir, "pediatric_all_patient_combined_clinical_input.rds")
pediatric_patient_clinical <- readRDS(pediatric_patient_clinical)

# transcriptomically similar samples table
fname <- file.path(output_dir, "transciptomically_similar_pediatric.rds")
transciptomically_similar_pediatric <- transciptomically_similar(all_cor = nn_table, 
                                                                 clin_data = pediatric_patient_clinical)

# save output
saveRDS(transciptomically_similar_pediatric, file = fname)