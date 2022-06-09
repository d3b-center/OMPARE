# define directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "coseq_detect")
output_dir <- file.path(patient_dir, "output", "coseq_detect", "adult")
dir.create(output_dir, showWarnings = F, recursive = T)

# run coseq_detect
coseq_script <- file.path(module_dir, "utils", 'coseq_detect.R')
cmd <- paste('Rscript', coseq_script,
             '--ref_cancer_dir', adult_cancer_dir, 
             '--output_dir', output_dir)
system(cmd)
