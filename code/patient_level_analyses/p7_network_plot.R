# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'network_plot.R'))

# network plot
network_plot_output <- network_plot(numGenes = 1000, 
                                    geneMania = gene_mania)

# save output
# saveRDS(network_plot_output, file = file.path(topDir, "output", "network_plot_output.rds"))