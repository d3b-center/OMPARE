
#### Define Directories --------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
module_dir <- file.path(root_dir, "code", "gsnca_analysis")
patient_dir <- file.path(root_dir, "results", patient)

input_dir <- file.path(patient_dir, "output", "transcriptomically_similar_analysis")
output_dir <- file.path(patient_dir, "output", "gsnca_analysis")
dir.create(output_dir, showWarnings = F, recursive = T)

# run coseq_detect
gsnca_script <- file.path(module_dir, 'gsnca_analysis.R')
cmd <- paste('Rscript', gsnca_script,
             '--patient', patient)
system(cmd)