
#### Define Directories --------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "gsnca_analysis")

# run gsnca
gsnca_script <- file.path(module_dir, 'gsnca_analysis.R')
cmd <- paste('Rscript', gsnca_script,
             '--patient', patient)
system(cmd)