# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'network_plot.R'))

# input transcriptomic drug recommendations output
fname <- file.path(topDir, "output", "transcriptome_drug_rec.rds")
transcriptome_drug_rec_output <- readRDS(fname)

# network plot
network_plot_output <- network_plot(transcriptome_drug_rec_output, geneMania = gene_mania)
