# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'diffexpr_genes_barplot.R'))

# input files
genes_up <- read.delim(file.path(topDir, "output", paste0(sampleInfo$subjectID, '_summary_DE_Genes_Up.txt')))
genes_down <- read.delim(file.path(topDir, "output", paste0(sampleInfo$subjectID, '_summary_DE_Genes_Down.txt')))

# call function
fname <- file.path(topDir, "output", "diffexpr_genes_barplot_output.rds")
diffexpr_genes_barplot_gtex <- diffexpr_genes_barplot(genes_up, genes_down, comparison_study = 'GTExBrain_1152', cancer_genes = cancer_genes$Gene_Symbol)
diffexpr_genes_barplot_pbta_hgg <- diffexpr_genes_barplot(genes_up, genes_down, comparison_study = 'PBTA_HGG_189', cancer_genes = cancer_genes$Gene_Symbol)
diffexpr_genes_barplot_pbta <- diffexpr_genes_barplot(genes_up, genes_down, comparison_study = 'PBTA_ALL_1035', cancer_genes = cancer_genes$Gene_Symbol)
diffexpr_genes_barplot_output <- list(diffexpr_genes_barplot_gtex = diffexpr_genes_barplot_gtex,
                                      diffexpr_genes_barplot_pbta_hgg = diffexpr_genes_barplot_pbta_hgg,
                                      diffexpr_genes_barplot_pbta = diffexpr_genes_barplot_pbta)
saveRDS(diffexpr_genes_barplot_output, file = fname)
