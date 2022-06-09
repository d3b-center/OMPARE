# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "p1_modules")
output_dir <- file.path(patient_dir, "output")
dir.create(output_dir, recursive = T, showWarnings = F)

# source functions
source(file.path(module_dir, "utils", 'genomic_summary.R'))

# call function
genomic_summary_output <- genomic_summary(key_clinical_findings_output = key_clinical_findings_output, 
                                          all_findings_output = all_findings_output,
                                          rnaseq_analysis_output = rnaseq_analysis_output)
saveRDS(genomic_summary_output, file = fname)