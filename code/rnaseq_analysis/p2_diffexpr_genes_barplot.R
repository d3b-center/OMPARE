# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "rnaseq_analysis")
output_dir <- file.path(patient_dir, "output", "rnaseq_analysis")
dir.create(output_dir, recursive = T, showWarnings = F)

# source functions
source(file.path(module_dir, "utils", 'diffexpr_genes_barplot.R'))

# cancer genes
cancer_genes <- readRDS(file.path(root_dir, "data", "cancer_gene_list.rds"))

# input files
genes_up <- read.delim(file.path(output_dir, paste0(patient, '_summary_DE_Genes_Up.txt')))
genes_down <- read.delim(file.path(output_dir, paste0(patient, '_summary_DE_Genes_Down.txt')))
genes_diffexpr <- rbind(genes_up, genes_down)

# call function
fname <- file.path(output_dir, "diffexpr_genes_barplot_output.pdf")
diffexpr_genes_barplot_output <- plyr::dlply(genes_diffexpr, 
                                             .variables = "comparison", 
                                             .fun = function(x) diffexpr_genes_barplot(x, cancer_genes = cancer_genes$Gene_Symbol))

# save to pdf
pdf(fname)
diffexpr_genes_barplot_output
dev.off()
