# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "transcriptomically_similar_analysis")
output_dir <- file.path(patient_dir, "output", "transcriptomically_similar_analysis")
dir.create(output_dir, showWarnings = F, recursive = T)

# source functions
source(file.path(module_dir, "utils", 'mutational_analysis.R'))
source(file.path(module_dir, "utils", 'recurrent_alterations_plots.R'))
source(file.path(module_dir, "utils", 'shared_alterations_plots.R'))

# load inputs
pediatric_cancer_dir <- file.path(data_dir, "pediatric_data")
nn_tpm_input <- file.path(output_dir, "pediatric_all_nn_tpm.rds")
nn_tpm_input <- readRDS(nn_tpm_input)
all_findings_output <- readRDS(file.path(patient_dir, "output", "all_findings_output.rds"))
key_clinical_findings_output <- readRDS(file.path(patient_dir, "output", "key_clinical_findings_output.rds"))
pediatric_patient_clinical <-  file.path(output_dir, "pediatric_all_patient_combined_clinical_input.rds")
pediatric_patient_clinical <- readRDS(pediatric_patient_clinical)

# recurrent alterations
fname <- file.path(output_dir, "mutational_analysis_pediatric.rds")
mutational_analysis_pediatric <- mutational_analysis(nn_tpm_input = nn_tpm_input, 
                                                     ref_cancer_dir = pediatric_cancer_dir,
                                                     all_findings_output = all_findings_output,
                                                     key_clinical_findings_output = key_clinical_findings_output,
                                                     clinical = pediatric_patient_clinical,
                                                     comparison = "pediatric")

# save output
saveRDS(mutational_analysis_pediatric, file = fname)

# convert to plots
# recurrent alterations
recurrent_alterations_plots(ref_cancer_dir = pediatric_cancer_dir, 
                            patient_maf = filtered_maf, 
                            patient_cnv = filtered_cnv,
                            mutational_analysis_output = mutational_analysis_pediatric, 
                            suffix = "pediatric")

# shared alterations
shared_alterations_plots(ref_cancer_dir = pediatric_cancer_dir, 
                         patient_maf = filtered_maf, 
                         patient_cnv = filtered_cnv,
                         mutational_analysis_output = mutational_analysis_pediatric, 
                         suffix = "pediatric")
