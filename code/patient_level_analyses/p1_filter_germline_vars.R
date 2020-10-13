# filtered germline variants

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'filter_germline_vars.R')) # filter germline variants

# filter germline data
filtered_germ_vars <- filter_germline_vars(mutData.germ = mutData.germ, myGermlineMarkers = germline_markers)

# save output
saveRDS(filtered_germ_vars, file = file.path(topDir, "output", "filtered_germline_vars.rds"))
