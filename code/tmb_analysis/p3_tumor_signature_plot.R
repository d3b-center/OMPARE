# Tumor signature plot
library(decompTumor2Sig)

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
module_dir <- file.path(root_dir, "code", "tmb_analysis")
output_dir <- file.path(patient_dir, "output", "tmb_analysis")
dir.create(output_dir, recursive = T, showWarnings = F)

# source functions
source(file.path(module_dir, "utils", 'tumor_signature_plot.R')) 

# output file
fname <- file.path(output_dir, "tumor_signature_output.rds")

# mutational signatures
signatures <- readAlexandrovSignatures(file.path(data_dir, 'signatures_probabilities.txt'))

# input maf file
patient_maf <- list.files(path = file.path(patient_dir, "simple-variants"), pattern = "consensus", full.names = T)
patient_maf <- data.table::fread(patient_maf)

# tumor signature plot
tumor_signature_output <- tumor_signature_plot(maf_data = patient_maf, 
                                               signatures = signatures, 
                                               output_dir = output_dir)

# save output
saveRDS(tumor_signature_output, file = fname)
