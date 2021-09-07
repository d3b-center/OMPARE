# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "rnaseq_analysis")
output_dir <- file.path(patient_dir, "output", "rnaseq_analysis")
dir.create(output_dir, recursive = T, showWarnings = F)

# source functions
source(file.path(module_dir, "utils", 'diffexpr_genes_barplot.R'))

# input files
genes_up <- read.delim(file.path(output_dir, paste0(patient, '_summary_DE_Genes_Up.txt')))
genes_down <- read.delim(file.path(output_dir, paste0(patient, '_summary_DE_Genes_Down.txt')))

# call function
fname <- file.path(output_dir, "diffexpr_genes_barplot_output.rds")
diffexpr_genes_barplot_gtex <- diffexpr_genes_barplot(genes_up, genes_down, comparison_study = 'GTExBrain_1152', cancer_genes = cancer_genes$Gene_Symbol)
diffexpr_genes_barplot_pbta_hgg <- diffexpr_genes_barplot(genes_up, genes_down, comparison_study = 'PBTA_HGG_189', cancer_genes = cancer_genes$Gene_Symbol)
diffexpr_genes_barplot_pbta <- diffexpr_genes_barplot(genes_up, genes_down, comparison_study = 'PBTA_ALL_1035', cancer_genes = cancer_genes$Gene_Symbol)
diffexpr_genes_barplot_output <- list(diffexpr_genes_barplot_gtex = diffexpr_genes_barplot_gtex,
                                      diffexpr_genes_barplot_pbta_hgg = diffexpr_genes_barplot_pbta_hgg,
                                      diffexpr_genes_barplot_pbta = diffexpr_genes_barplot_pbta)
saveRDS(diffexpr_genes_barplot_output, file = fname)
