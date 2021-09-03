# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
module_dir <- file.path(root_dir, "code", "genomic_landscape_plots")
output_dir <- file.path(patient_dir, "output", "genomic_landscape_plots")
dir.create(output_dir, showWarnings = F, recursive = T)

# source functions
source(file.path(module_dir, "utils", "network_plot.R"))

# input transcriptomic drug recommendations output
fname <- file.path(root_dir, "output", "drug_recommendations", "transcriptome_drug_rec.rds")
transcriptome_drug_rec_output <- readRDS(fname)

# gene mania reference
gene_mania <- read.delim(file.path(data_dir, 'GeneManiaNetwork.txt'), stringsAsFactors =F)

# network plot
network_plot_output <- network_plot(transcriptome_drug_rec_output, 
                                    geneMania = gene_mania,
                                    filtered_mutations = mutDataFilt, 
                                    filtered_fusions = fusData, 
                                    rnaseq_analysis_output = rnaseq_analysis_output)
