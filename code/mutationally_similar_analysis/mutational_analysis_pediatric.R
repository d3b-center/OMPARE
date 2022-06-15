# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "mutationally_similar_analysis")
output_dir <- file.path(patient_dir, "output", "mutationally_similar_analysis")
dir.create(output_dir, showWarnings = F, recursive = T)

# source functions
source(file.path(module_dir, "utils", 'mutational_analysis.R'))
source(file.path(module_dir, "utils", 'recurrent_alterations_plots.R'))
source(file.path(module_dir, "utils", 'shared_alterations_plots.R'))

# load inputs
nn_table <- file.path(patient_dir, "output", "mut_distance_calc", "pediatric", "nn_table.rds")
nn_table <- readRDS(nn_table)
sample_id_interest <- sample_info %>%
  pull(sample_id)

# recurrent alterations
fname <- file.path(output_dir, "mutational_analysis_pediatric.rds")
mutational_analysis_pediatric <- mutational_analysis(nn_table = nn_table, 
                                                     sample_id_interest = sample_id_interest,
                                                     ref_cancer_dir = pediatric_cancer_dir,
                                                     all_findings_output = all_findings_output,
                                                     key_clinical_findings_output = key_clinical_findings_output,
                                                     comparison = "pediatric")

# save output
saveRDS(mutational_analysis_pediatric, file = fname)

# convert to plots
# recurrent alterations
recurrent_alterations_plots(ref_cancer_dir = pediatric_cancer_dir, 
                            patient_maf = filtered_maf, 
                            patient_cnv = NULL,
                            mutational_analysis_output = mutational_analysis_pediatric, 
                            prefix = "pediatric")

# shared alterations
shared_alterations_plots(ref_cancer_dir = pediatric_cancer_dir, 
                         patient_maf = filtered_maf, 
                         patient_cnv = NULL,
                         mutational_analysis_output = mutational_analysis_pediatric, 
                         prefix = "pediatric")
