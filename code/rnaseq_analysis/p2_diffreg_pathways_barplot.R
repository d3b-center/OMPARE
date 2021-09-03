# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "rnaseq_analysis")
output_dir <- file.path(patient_dir, "output", "rnaseq_analysis")
dir.create(output_dir, recursive = T, showWarnings = F)

# source functions
source(file.path(module_dir, "utils", 'diffreg_pathways_barplot.R'))

# input files
pathways_up <- read.delim(file.path(output_dir, paste0(patient, '_summary_Pathways_Up.txt')))
pathways_down <- read.delim(file.path(output_dir, paste0(patient, '_summary_Pathways_Down.txt')))

# call function
fname <- file.path(output_dir, "diffreg_pathways_barplot_output.rds")
if(!file.exists(fname)){
  diffreg_pathways_barplot_gtex <- diffreg_pathways_barplot(pathways_up, pathways_down, comparison_study = 'GTExBrain_1152')
  diffreg_pathways_barplot_pbta_hgg <- diffreg_pathways_barplot(pathways_up, pathways_down, comparison_study = 'PBTA_HGG_189')
  diffreg_pathways_barplot_pbta <- diffreg_pathways_barplot(pathways_up, pathways_down, comparison_study = 'PBTA_ALL_1035')
  diffreg_pathways_barplot_output <- list(diffreg_pathways_barplot_gtex = diffreg_pathways_barplot_gtex,
                                          diffreg_pathways_barplot_pbta_hgg = diffreg_pathways_barplot_pbta_hgg,
                                          diffreg_pathways_barplot_pbta = diffreg_pathways_barplot_pbta)
  saveRDS(diffreg_pathways_barplot_output, file = fname)
} else {
  diffreg_pathways_barplot_output <- readRDS(fname)
}
