# Author: Komal S. Rathi
# Function to define directories

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
fig_dir <- file.path(root_dir, 'figures')
code_dir <- file.path(root_dir, 'code')
utils_dir <- file.path(code_dir, 'utils')
data_dir <- file.path(root_dir, 'data')
ref_dir <- file.path(data_dir, 'Reference')
ref_input_file_dir <- file.path(ref_dir, "input_files")
