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
genes_up <- read.delim(file.path(output_dir, 'genes_up.txt'))
genes_down <- read.delim(file.path(output_dir, 'genes_down.txt'))

# call function
comparisons <- unique(genes_up$comparison)
diffexpr_genes_barplot_output <- list()
for(i in 1:length(comparisons)){
  diffexpr_genes_barplot_output[[i]] <- diffexpr_genes_barplot(genes_up, genes_down, 
                                                               comparison_study = comparisons[i], 
                                                               cancer_genes = cancer_genes$Gene_Symbol)
}
fname <- file.path(output_dir, "diffexpr_genes_barplot_output.pdf")
pdf(file = fname, width = 14, height = 6)
p <- ggpubr::ggarrange(plotlist = diffexpr_genes_barplot_output, nrow = 1)
print(p)
dev.off()
