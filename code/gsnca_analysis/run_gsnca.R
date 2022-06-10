
# define directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "gsnca_analysis")
output_dir <- file.path(patient_dir, "output", "gsnca_analysis")
dir.create(output_dir, showWarnings = F, recursive = T)

# sample_info is sample info on patient of interest
patient_of_interest <- sample_info %>%
  filter(experimental_strategy == "RNA-Seq") %>%
  pull(Kids_First_Biospecimen_ID)

# coseq clustered data (pediatric tumors)
nb_cluster_assigned <- file.path(patient_dir, "output", "coseq_detect", "pediatric", "cancer_group_of_interest_nb_cluster_assigned.tsv")

# directories
pediatric_cancer_dir <- file.path(root_dir, "data", "pediatric_data")
normal_tissue_dir <- file.path(root_dir, "data", "normal_data")
adult_cancer_dir <- file.path(data_dir, "adult_data")

# run gsnca
gsnca_script <- file.path(module_dir, 'gsnca_analysis.R')
cmd <- paste('Rscript', gsnca_script,
             '--patient', patient_of_interest,
             '--pediatric_cancer_dir', pediatric_cancer_dir, 
             '--norm_data_dir', normal_tissue_dir, 
             '--adult_data_dir', adult_cancer_dir,
             '--nb_cluster_assigned', nb_cluster_assigned,
             '--output_dir', output_dir)
system(cmd)
