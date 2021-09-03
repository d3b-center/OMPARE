# key clinical findings

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "p1_modules")

# source functions
source(file.path(module_dir, "utils", 'key_clinical_findings.R'))

# call function
key_clinical_findings_output <- key_clinical_findings(snv_caller = snv_caller, 
                                                      all_findings_output = all_findings_output)
