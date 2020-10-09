# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'cnv_plot.R'))

# cnv plot
cnv_plot(myCnvData = cnvRatioData, fname = file.path(topDir, "output", "cnv_plot.png"))
