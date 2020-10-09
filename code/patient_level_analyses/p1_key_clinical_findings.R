# key clinical findings

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'key_clinical_findings.R'))

# call function
key_clinical_findings_output <- key_clinical_findings(snv_pattern = snv_pattern, 
                                                      all_findings_output = all_findings_output)
