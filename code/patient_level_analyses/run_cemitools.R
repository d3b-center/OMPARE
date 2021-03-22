# library(tidyverse)
library(sva)
library(CEMiTool)
library(tidyverse)
library(dplyr)

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))
source(file.path(patient_level_analyses_utils, "data_formatting_functions.R"))

# depends on this version of hgg-dmg data
hgg_dmg_20201202 <- file.path(ref_dir, 'hgg-dmg-integration', '20201202-data')
pnoc008_dir <- file.path(ref_dir, 'pnoc008')
patient <- sampleInfo$subjectID

# output directory
patient_output_dir <- file.path(results_dir, patient, "output")
cemitools_dir <- file.path(results_dir, patient, "CEMITools")
dir.create(cemitools_dir, showWarnings = F, recursive = T)

# ge based clustering
ge_based_clustering <- file.path(hgg_dmg_20201202, 'CC_based_heatmap_Distance_euclidean_finalLinkage_average_clusterAlg_KM_expct_counts_VST_cluster_and_annotation.tsv')
ge_based_clustering <- read.delim(ge_based_clustering)
ge_based_clustering <- ge_based_clustering %>%
  dplyr::select(Sample.Names, CC)

# pbta histology
pbta_histology <- read.delim(file.path(hgg_dmg_20201202, 'pbta-histologies.tsv'))
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
  tibble::remove_rownames() %>%
  tibble::column_to_rownames("Kids_First_Biospecimen_ID")
  
# hgg-dmg expected counts
exp_counts <- readRDS(file.path(hgg_dmg_20201202, 'pbta-hgat-dx-prog-pm-gene-counts-rsem-expected_count-uncorrected.rds'))
exp_counts <- exp_counts[rownames(pnoc008_counts),]

# merge both and subset to 122 ge based clustered samples + pnoc008 patient of interest
counts_mat <- cbind(pnoc008_counts, exp_counts)
counts_mat <- counts_mat[,rownames(clinical)]

# batch correct using Combat_seq (used for counts)
exp_counts_corrected <- ComBat_seq(counts = as.matrix(counts_mat), batch = clinical$batch)
saveRDS(exp_counts_corrected, file = file.path(cemitools_dir, "expected_counts_corrected.rds"))

# umap and get top 20 nearest neighbors
exp_most_var <- get_most_variable_for_umap(expr_corrected = exp_counts_corrected)
umap_output <- get_umap_output(expr_most_var = exp_most_var)
sampleInfo <- data.frame(subjectID = pnoc008_clinical$Kids_First_Biospecimen_ID)
nn_table <- extract_umap_nearest_neighbor_table(umap_out = umap_output, expr_most_var = exp_most_var, sampleInfo = sampleInfo)

# get clusters corresponding to nearest neighbors and assign max occurence to pnoc008 patient of interest
x <- clinical[nn_table$nearest_neighbor,'CC']
pnoc008_cluster <- as.numeric(names(table(x))[which.max(table(x))])
clinical[pnoc008_clinical$Kids_First_Biospecimen_ID,'CC'] <- pnoc008_cluster

# subset to Kids_First_Biospecimen_ID and CC column for CEMItools
clinical <- clinical %>%
  tibble::rownames_to_column('Kids_First_Biospecimen_ID') %>%
  dplyr::select(Kids_First_Biospecimen_ID, CC)
saveRDS(clinical, file = file.path(cemitools_dir, "clustered_samples.rds"))

# run CEMItools
n <- 100
cem <- cemitool(as.data.frame(exp_counts_corrected), 
                clinical, 
                filter = T,
                cor_functio = 'bicor', 
                network_type = 'signed',
                tom_type = 'signed',
                sample_name_column = 'Kids_First_Biospecimen_ID',
                class_column = 'CC',
                merge_similar = T,
                apply_vst = T)

# write out hubs and summary
hubs <- get_hubs(cem,n)
summary <- mod_summary(cem)
saveRDS(hubs, file = file.path(cemitools_dir, 'hubs.rds'))
saveRDS(summary, file = file.path(cemitools_dir, 'summary.rds'))

# generate heatmap of gene set enrichment analysis
cem <- mod_gsea(cem)
cem <- plot_gsea(cem)

# plot gene expression within each module
cem <- plot_profile(cem)

# read GMT file - reactome file
gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_in <- read_gmt(gmt_fname)

# perform over representation analysis
cem <- mod_ora(cem, gmt_in)

# plot ora results
cem <- plot_ora(cem)

# read interactions
int_fname <- system.file("extdata", "interactions.tsv", package = "CEMiTool")
int_df <- read.delim(int_fname)

# plot interactions
interactions_data(cem) <- int_df # add interactions
cem <- plot_interactions(cem) # generate plot

# save interaction plots for all modules
r(function(x, y) { CEMiTool::diagnostic_report(cem = x, directory = y, force = T) }, args = list(cem, cemitools_dir))
# diagnostic_report(cem, directory = cemitools_dir, force = T)

# output
r(function(x, y) { CEMiTool::generate_report(cem = x, directory = y, force = T) }, args = list(cem, cemitools_dir))
# generate_report(cem, directory = cemitools_dir, force = T)
write_files(cem, directory = cemitools_dir, force = T)
save_plots(cem, "all", directory = cemitools_dir, force = T)

# positively correlated modules for pnoc008 cluster
corr_modules <- cem@enrichment$es
corr_modules <- corr_modules %>%
  filter(get(as.character(pnoc008_cluster)) > 0) %>%
  .$pathway

# genes-module mapping for positively correlated modules
corr_modules <- module_genes(cem) %>%
  filter(modules %in% corr_modules)

# annotate targetable hubs
fname <- file.path(patient_output_dir, "transcriptome_drug_rec.rds")
dge_genes <- readRDS(fname)
dge_genes <- dge_genes %>% 
  mutate(Hub_Gene = ifelse(Comparison == "PBTA_HGG_182" & Gene %in% corr_modules$genes, TRUE, FALSE))

# rewrite transcriptiomic based drug recommendations
saveRDS(dge_genes, file = fname)
