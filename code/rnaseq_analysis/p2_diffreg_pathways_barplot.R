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
pathways_diffreg <- rbind(pathways_up, pathways_down)

# call function
fname <- file.path(output_dir, "diffreg_pathways_barplot_output.pdf")
diffreg_pathways_barplot_output <- plyr::dlply(pathways_diffreg, 
                                             .variables = "comparison", 
                                             .fun = function(x) diffreg_pathways_barplot(x))

# save to pdf
pdf(fname, width = 12)
diffreg_pathways_barplot_output
dev.off()

