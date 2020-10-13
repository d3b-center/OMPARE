# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'disease_specific_information.R'))

# call function
disease_specific_information_output <- disease_specific_information(disease_specific_fields = disease_specific_fields,
                                                                    snv_pattern = snv_pattern,  
                                                                    all_findings_output = all_findings_output)