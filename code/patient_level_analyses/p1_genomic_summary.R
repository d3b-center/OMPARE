# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'genomic_summary.R'))

# call function
genomic_summary_output <- genomic_summary(snv_caller = snv_caller, 
                                          key_clinical_findings_output = key_clinical_findings_output, 
                                          all_findings_output = all_findings_output)
