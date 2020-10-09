# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# prepare input files for oncogrid
source(file.path(patient_level_analyses_utils, 'prepare_files_oncogrid.R'))

# plot oncogrid
source(file.path(patient_level_analyses_utils, 'plot_oncogrid.R'))
