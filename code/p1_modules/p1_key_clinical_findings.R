# key clinical findings

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "p1_modules")
output_dir <- file.path(patient_dir, "output")
dir.create(output_dir, recursive = T, showWarnings = F)

# input data
if(file.exists(file.path(output_dir, "rnaseq_analysis", "rnaseq_analysis_output.rds"))){
  rnaseq_analysis_output <- readRDS(file.path(output_dir, "rnaseq_analysis", "rnaseq_analysis_output.rds"))
} else {
  rnaseq_analysis_output <- NULL
}
all_findings_output <- readRDS(file.path(output_dir, "all_findings_output.rds"))

# source functions
source(file.path(module_dir, "utils", 'key_clinical_findings.R'))

# call function
fname <- file.path(output_dir, "key_clinical_findings_output.rds")
key_clinical_findings_output <- key_clinical_findings(rnaseq_analysis_output = rnaseq_analysis_output,
                                                      all_findings_output = all_findings_output)
saveRDS(key_clinical_findings_output, file = fname)
