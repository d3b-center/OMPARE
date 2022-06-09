# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "oncogrid_analysis")
output_dir <- file.path(patient_dir, "output", "oncogrid_analysis")
dir.create(output_dir, showWarnings = F, recursive = T)

# source functions
source(file.path(module_dir, "utils", "plot_oncogrid.R"))

# prepare input files for oncogrid
source(file.path(module_dir, "utils", "prepare_files_oncogrid.R"))
prepare_files_oncogrid(pediatric_cancer_dir = pediatric_cancer_dir,
                       sample_info = sample_info,
                       output_dir = output_dir)

# plot oncogrid
fname <- file.path(output_dir, "complexheatmap_oncogrid.pdf")
plot_oncogrid(input_dir = output_dir, output_file = fname, width = 40, height = 14)
