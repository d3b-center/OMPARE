# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "transcriptomically_similar_analysis")
output_dir <- file.path(patient_dir, "output", "transcriptomically_similar_analysis")
dir.create(output_dir, showWarnings = F, recursive = T)

# source functions
source(file.path(module_dir, "utils", 'ssgsea.R'))
source(file.path(module_dir, "utils", 'plot_ssgsea.R'))

# compute ssgsea 
fname <- file.path(output_dir, "ssgsea_scores_pediatric.rds")
ssgsea_pediatric <- ssgsea(top_cor = pbta_pnoc008_nn_tpm, patient_of_interest = patient)
saveRDS(ssgsea_pediatric, file = fname)

# plot ssgea
fname <- file.path(output_dir, 'ssgsea_scores_pediatric.pdf')
ssgsea_pediatric <- plot_ssgsea(ssgsea_pediatric)
ggsave(filename = fname, plot = ssgsea_pediatric, 
       device = "pdf", width = 15, height = 14)
