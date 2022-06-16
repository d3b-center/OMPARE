# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
module_dir <- file.path(root_dir, "code", "genomic_landscape_plots")
output_dir <- file.path(patient_dir, "output", "genomic_landscape_plots")
dir.create(output_dir, showWarnings = F, recursive = T)

# source functions
source(file.path(module_dir, "utils", "circos_plot.R"))

# cancer Genes (annoFuse)
cancer_genes <- readRDS(file.path(data_dir, 'cancer_gene_list.rds'))

# chr coordinates to gene symbol map
chr_map <- read.delim(file.path(data_dir, 'mart_export_genechr_mapping.txt'), stringsAsFactors = F, check.names = F)
colnames(chr_map) <- c("hgnc_symbol", "gene_start", "gene_end", "chromosome")

# read data
rnaseq_analysis_output <- readRDS(file.path(patient_dir, "output", "rnaseq_analysis", "rnaseq_analysis_output.rds"))

# circos plot
fname <- file.path(file.path(output_dir, "circos_plot.png"))
circos_plot(chr_map = chr_map, 
            cancer_genes = cancer_genes, 
            fname = fname,
            filtered_mutations = filtered_maf, 
            rnaseq_analysis_output = rnaseq_analysis_output, 
            filtered_fusions = filtered_fusions)
