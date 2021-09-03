# filter germline variants

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "p1_modules")
output_dir <- file.path(patient_dir, "output")
dir.create(output_dir, recursive = T, showWarnings = F)

# source functions
source(file.path(module_dir, "utils", 'filter_germline_vars.R')) # filter germline variants

# this is snv caller independent so run only if output does not exist
fname <- file.path(output_dir, "filtered_germline_vars.rds")
if(!file.exists(fname)){
  # Germline markers
  germline_markers <- readRDS(file.path(data_dir, 'germline_markers_list.rds'))
  
  # filter germline data
  filtered_germ_vars <- filter_germline_vars(patient_dir, germline_markers = germline_markers)
  
  # save output
  saveRDS(filtered_germ_vars, file = fname)
} else {
  filtered_germ_vars <- readRDS(fname)
}
