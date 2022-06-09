# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "rnaseq_analysis")
output_dir <- file.path(patient_dir, "output", "rnaseq_analysis")
dir.create(output_dir, recursive = T, showWarnings = F)

# source functions
source(file.path(module_dir, "utils", 'diffreg_pathways_barplot.R'))

# input files
pathways_up <- read.delim(file.path(output_dir, 'pathways_up.txt'))
pathways_down <- read.delim(file.path(output_dir, 'pathways_down.txt'))

# call function
fname <- file.path(output_dir, "diffreg_pathways_barplot_output.pdf")

comparisons <- unique(pathways_up$comparison)
diffreg_pathways_barplot_output <- list()
pdf(file = fname, width = 14, height = 6)
for(i in 1:length(comparisons)){
  p <- diffreg_pathways_barplot(pathways_up, pathways_down, 
                                comparison_study = comparisons[i])
  print(p)
}
dev.off()
