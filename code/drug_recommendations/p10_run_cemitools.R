suppressPackageStartupMessages({
  library(sva)
  library(tidyverse)
  library(dplyr)
  library(optparse)
})

# arguments
option_list <- list(
  make_option(c("--patient"), type = "character",
              help = "Patient identifier (PNOC008-22, C3342894...)")
)
opt <- parse_args(OptionParser(option_list = option_list))
patient <- opt$patient

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
module_dir <- file.path(root_dir, "code", "drug_recommendations")

# output directory
output_dir <- file.path(root_dir, "results", patient, "output", "drug_recommendations")
cemitools_dir <- file.path(output_dir, "CEMITools")
dir.create(cemitools_dir, showWarnings = F, recursive = T)

# source functions
source(file.path(root_dir, "code", "transcriptomically_similar_analysis", "utils", "get_most_variable_for_umap.R"))
source(file.path(root_dir, "code", "transcriptomically_similar_analysis", "utils", "get_umap_output.R"))
source(file.path(root_dir, "code", "transcriptomically_similar_analysis", "utils", "extract_umap_nearest_neighbor_table.R"))
source(file.path(root_dir, "code", "utils", "pubTheme.R"))
source(file.path(module_dir, "utils", "run_cemitools.R"))

# depends on 20201202 version of hgg-dmg data
hgg_dmg_20201202 <- file.path(data_dir, "hgg-dmg-integration", "20201202-data")
pnoc008_dir <- file.path(data_dir, "pnoc008")

# ge based clustering
ge_based_clustering <- file.path(hgg_dmg_20201202, "CC_based_heatmap_Distance_euclidean_finalLinkage_average_clusterAlg_KM_expct_counts_VST_cluster_and_annotation.tsv")
ge_based_clustering <- read.delim(ge_based_clustering)
ge_based_clustering <- ge_based_clustering %>%
  dplyr::select(Sample.Names, CC)

# pbta histology
pbta_histology <- read.delim(file.path(hgg_dmg_20201202, "pbta-histologies.tsv"))
pbta_histology <- pbta_histology %>%
  filter(experimental_strategy == "RNA-Seq",
         pbta.hgat.dx == "Y",
         cohort_participant_id %in% ge_based_clustering$Sample.Names) %>%
  mutate(batch = paste0(cohort, "_", RNA_library)) %>%
  dplyr::select(Kids_First_Biospecimen_ID, cohort_participant_id, batch)

# pnoc008 clinical file (subset to patient of interest)
pnoc008_clinical <- readRDS(file.path(pnoc008_dir, 'pnoc008_clinical.rds'))
pnoc008_clinical <- pnoc008_clinical %>%
  filter(subjectID == patient) %>%
  mutate(batch = paste0(study_id, "_", library_name)) %>%
  dplyr::select(Kids_First_Biospecimen_ID, cohort_participant_id, batch)

# pnoc008 counts (subset to patient of interest)
pnoc008_counts <- readRDS(file.path(pnoc008_dir, 'pnoc008_counts_matrix.rds'))
pnoc008_counts <- pnoc008_counts[,patient, drop = F]
colnames(pnoc008_counts) <- pnoc008_clinical$Kids_First_Biospecimen_ID

# merge clinical file
clinical <- pnoc008_clinical %>%
  rbind(pbta_histology) %>%
  left_join(ge_based_clustering, by = c("cohort_participant_id" = "Sample.Names")) %>%
  unique() %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames("Kids_First_Biospecimen_ID")

# hgg-dmg expected counts
exp_counts <- readRDS(file.path(hgg_dmg_20201202, "pbta-hgat-dx-prog-pm-gene-counts-rsem-expected_count-uncorrected.rds"))
rows_to_keep <- intersect(rownames(exp_counts), rownames(pnoc008_counts))
exp_counts <- exp_counts[rows_to_keep,]
pnoc008_counts <- pnoc008_counts[rows_to_keep,,drop = F]

# merge both and subset to 122 ge-based clustered samples + pnoc008 patient of interest
counts_mat <- cbind(pnoc008_counts, exp_counts)
counts_mat <- counts_mat[,rownames(clinical)]

# batch correct using Combat_seq (used for counts)
if(min(table(clinical$batch)) == 1){
  exp_counts_corrected <- ComBat(dat = as.matrix(log2(counts_mat + 1)), batch = clinical$batch)
  exp_counts_corrected <- 2^exp_counts_corrected
} else {
  exp_counts_corrected <- ComBat_seq(counts = as.matrix(counts_mat), batch = clinical$batch)
}
saveRDS(exp_counts_corrected, file = file.path(cemitools_dir, "expected_counts_corrected.rds"))

# umap and get top 20 nearest neighbors
exp_most_var <- get_most_variable_for_umap(expr_corrected = exp_counts_corrected)
umap_output <- get_umap_output(expr_most_var = exp_most_var)
sample_info <- data.frame(subjectID = pnoc008_clinical$Kids_First_Biospecimen_ID)
nn_table <- extract_umap_nearest_neighbor_table(umap_out = umap_output, 
                                                expr_most_var = exp_most_var, 
                                                patient_of_interest = sample_info$subjectID)
saveRDS(umap_output, file = file.path(cemitools_dir, "umap_output.rds"))
saveRDS(nn_table, file = file.path(cemitools_dir, "umap_top_20_neighbors_output.rds"))

# get clusters corresponding to nearest neighbors and assign max occurence to pnoc008 patient of interest
x <- clinical[nn_table$nearest_neighbor,'CC']
pnoc008_cluster <- as.numeric(names(table(x))[which.max(table(x))])
clinical[pnoc008_clinical$Kids_First_Biospecimen_ID,'CC'] <- pnoc008_cluster

# subset to Kids_First_Biospecimen_ID and CC column for CEMItools
clinical <- clinical %>%
  tibble::rownames_to_column('Kids_First_Biospecimen_ID') %>%
  dplyr::select(Kids_First_Biospecimen_ID, CC)
saveRDS(clinical, file = file.path(cemitools_dir, "clustered_samples.rds"))

# call function
run_cemitools(exp_counts_corrected = exp_counts_corrected, 
              clinical = clinical, 
              output_dir = output_dir, 
              cemitools_dir = cemitools_dir, 
              pnoc008_cluster = pnoc008_cluster)
