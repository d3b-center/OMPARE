# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'circos_plot.R'))

# circos plot
circos_plot(topDir = topDir, chrMap = chrMap, fname = file.path(topDir, "output", "circos_plot.png"))
