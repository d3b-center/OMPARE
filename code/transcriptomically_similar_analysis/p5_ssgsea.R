# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "transcriptomically_similar_analysis")
output_dir <- file.path(patient_dir, "output", "transcriptomically_similar_analysis")
dir.create(output_dir, showWarnings = F, recursive = T)

# source functions
source(file.path(module_dir, "utils", 'ssgsea.R'))
source(file.path(module_dir, "utils", 'plot_ssgsea.R'))

# load inputs
nn_tpm_input <- file.path(patient_dir, "output", "transcriptomically_similar_analysis", "pediatric_all_nn_tpm.rds")
nn_tpm_input <- readRDS(nn_tpm_input)

# compute ssgsea 
fname <- file.path(output_dir, "ssgsea_scores_pediatric.rds")
patient_of_interest <- sample_info %>%
  filter(experimental_strategy == "RNA-Seq") %>%
  .$Kids_First_Biospecimen_ID
ssgsea_pediatric <- ssgsea(nn_tpm_input = nn_tpm_input, patient_of_interest = patient_of_interest)
saveRDS(ssgsea_pediatric, file = fname)

# plot ssgea
fname <- file.path(output_dir, 'ssgsea_scores_pediatric.pdf')
ssgsea_pediatric_plot <- plot_ssgsea(ssgsea_pediatric)
ggsave(filename = fname, plot = ssgsea_pediatric_plot, 
       device = "pdf", width = 15, height = 14)
