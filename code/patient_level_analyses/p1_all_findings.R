# all findings table

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'all_findings.R'))

# output
all_findings_output <- all_findings(snv_caller = snv_caller)

