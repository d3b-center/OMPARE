# adult immune profile

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "immune_analysis")
output_dir <- file.path(patient_dir, "output", "immune_analysis")
dir.create(output_dir, recursive = T, showWarnings = F)

# source functions
source(file.path(module_dir, "utils", 'immune_profile.R'))
source(file.path(module_dir, "utils", 'plot_immune_profile.R'))

# load inputs
adult_immune_profile_input <- file.path(patient_dir, "output", "transcriptomically_similar_analysis", "adult_immune_profile_input.rds")
adult_immune_profile_input <- readRDS(adult_immune_profile_input)

# immune profile and save scores as well
fname <- file.path(output_dir, 'immune_scores_adult.rds')
adult_immune_profile <- immune_profile(fullmat = adult_immune_profile_input)
saveRDS(adult_immune_profile, file = fname)

# plot immune scores
fname <- file.path(output_dir, 'immune_scores_adult.pdf')
patient_of_interest <- sample_info %>%
  filter(experimental_strategy  == "RNA-Seq") %>%
  .$Kids_First_Biospecimen_ID
adult_immune_profile <- plot_immune_profile(xcell_scores = adult_immune_profile, patient_of_interest = patient_of_interest)
ggsave(filename = fname, plot = adult_immune_profile, 
       device = "pdf", width = 8, height = 10)
