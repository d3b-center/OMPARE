# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "drug_recommendations")

# output directory
output_dir <- file.path(patient_dir, "output", "drug_recommendations")
cemitools_dir <- file.path(output_dir, "CEMITools")
dir.create(cemitools_dir, showWarnings = F, recursive = T)

# source functions
source(file.path(root_dir, "code", "utils", "pubTheme.R"))
source(file.path(module_dir, "utils", "run_cemitools.R"))

# input files
coseq_output <- file.path(patient_dir, "output", "coseq_detect", "pediatric", "cancer_group_of_interest_nb_cluster_assigned.tsv")
coseq_output <- readr::read_tsv(coseq_output)
expr_mat <- readRDS(file.path(patient_dir, "output", "coseq_detect", "pediatric", "expression_matrix.rds"))
patient_of_interest <- sample_info %>%
  filter(experimental_strategy == "RNA-Seq") %>%
  pull(Kids_First_Biospecimen_ID)
poi_cluster <- coseq_output %>% 
  filter(Kids_First_Biospecimen_ID == patient_of_interest) %>% 
  pull(cluster_assigned_nb)
transcriptome_drug_rec_output <- readRDS(file.path(output_dir, "transcriptome_drug_rec.rds"))

# call function
run_cemitools(expr_mat = expr_mat, 
              annot = coseq_output, 
              output_dir = output_dir, 
              cemitools_dir = cemitools_dir, 
              poi_cluster = poi_cluster,
              transcriptome_drug_rec_output = transcriptome_drug_rec_output)
