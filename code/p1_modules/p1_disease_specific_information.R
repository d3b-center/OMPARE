# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "p1_modules")
output_dir <- file.path(patient_dir, "output")
dir.create(output_dir, recursive = T, showWarnings = F)

# source functions
source(file.path(module_dir, "utils", 'disease_specific_information.R'))

# HGG-specific genes
disease_specific_fields <- read.delim(file.path(root_dir, "data", 'DiseaseSpecificFields.txt'))

# call function
fname <- file.path(output_dir, "disease_specific_information_output.rds")
disease_specific_information_output <- disease_specific_information(disease_specific_fields = disease_specific_fields,
                                                                    snv_caller = snv_caller,  
                                                                    all_findings_output = all_findings_output)
saveRDS(disease_specific_information_output, file = fname)