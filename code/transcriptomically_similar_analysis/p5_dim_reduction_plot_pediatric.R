# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "transcriptomically_similar_analysis")
output_dir <- file.path(patient_dir, "output", "transcriptomically_similar_analysis")
dir.create(output_dir, showWarnings = F, recursive = T)

# source functions
source(file.path(module_dir, "utils", "dim_reduction_plot.R"))

# load inputs
umap_embedding <- file.path(output_dir, "pediatric_all_umap_embedding.rds")
umap_embedding <- readRDS(umap_embedding)
pediatric_patient_clinical <-  file.path(output_dir, "pediatric_all_patient_combined_clinical_input.rds")
pediatric_patient_clinical <- readRDS(pediatric_patient_clinical)
patient_of_interest <- sample_info %>%
  filter(experimental_strategy == "RNA-Seq") %>%
  .$Kids_First_Biospecimen_ID

# umap clustering plot
fname <- file.path(output_dir, "dim_reduction_plot_pediatric.rds")
dim_reduction_plot_pediatric <- dim_reduction_plot(umap_embedding = umap_embedding,
                                                   clindata = pediatric_patient_clinical,
                                                   patient_of_interest = patient_of_interest)
# save output
saveRDS(dim_reduction_plot_pediatric, file = fname)
