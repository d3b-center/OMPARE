# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "p1_modules")

# source functions
source(file.path(module_dir, "utils", 'genomic_summary.R'))

# call function
genomic_summary_output <- genomic_summary(snv_caller = snv_caller, 
                                          key_clinical_findings_output = key_clinical_findings_output, 
                                          all_findings_output = all_findings_output)
