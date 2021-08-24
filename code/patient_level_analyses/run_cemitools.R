# library(tidyverse)
library(sva)
library(CEMiTool)
library(tidyverse)
library(patchwork)
library(dplyr)

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))
source(file.path(patient_level_analyses_utils, "data_formatting_functions.R"))
source(file.path(patient_level_analyses_utils, "pubTheme.R"))

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
  unique() %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames("Kids_First_Biospecimen_ID")
  
# hgg-dmg expected counts
exp_counts <- readRDS(file.path(hgg_dmg_20201202, 'pbta-hgat-dx-prog-pm-gene-counts-rsem-expected_count-uncorrected.rds'))
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
nn_table <- extract_umap_nearest_neighbor_table(umap_out = umap_output, expr_most_var = exp_most_var, sampleInfo = sample_info)
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
                apply_vst = T, 
                verbose = F)

# write out hubs and summary
hubs <- get_hubs(cem, n, method = "kME")
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
gmt_in <- CEMiTool::read_gmt(gmt_fname)

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

# get pos/neg correlated modules for pnoc008 cluster using p-adj < 0.05 cut-off
corr_modules <- cem@enrichment$padj %>%
  filter(get(as.character(pnoc008_cluster)) < 0.05) %>%
    .$pathway
corr_modules <- cem@enrichment$nes %>%
  filter(pathway %in% corr_modules,
         pathway != "Not.Correlated") %>%
  mutate(direction = ifelse(get(as.character(pnoc008_cluster)) > 0, "pos", "neg"))

# get hub genes for pos/neg correlated modules with a cutoff of 0.5
# note: when using method = "kME", the hubs object does not behave like a normal list and so I am unable to use stack to unlist the list recursively
network_hubs <- stack(unlist(hubs, use.names = T))
network_hubs$ind <- gsub('Not.Correlated','Not_Correlated', network_hubs$ind)
network_hubs <-  cbind(values = network_hubs$values, reshape2::colsplit(network_hubs$ind, pattern = '\\.', names = c("module", "genes")))
network_hubs <- network_hubs %>%
  filter(module %in% corr_modules$pathway,
         values >= 0.5)
  
# annotate targetable hubs
fname <- file.path(patient_output_dir, "transcriptome_drug_rec.rds")
dge_genes <- readRDS(fname)
dge_genes <- dge_genes %>% 
  mutate(Network_Hub = ifelse(Comparison == "PBTA_HGG_189" & Gene %in% network_hubs$genes, "Yes", "No"))

# rewrite transcriptiomic based drug recommendations
saveRDS(dge_genes, file = fname)

# get ora data for pos/neg correlated modules
ora_dat <- ora_data(cem = cem)
ora_dat <- ora_dat %>%
  inner_join(corr_modules %>% dplyr::select(pathway, direction), by = c("Module" = "pathway"))
ora_dat <- ora_dat %>%
  group_by(Module) %>%
  mutate(order_val = row_number()) %>%
  arrange(p.adjust) %>%
  slice_head(n = 10) 
modules <- unique(ora_dat$Module)
p <- list()
for(i in 1:length(modules)){
  tmp <- ora_dat %>%
    filter(Module %in% modules[i])
  title <- paste("Module: ", unique(tmp$Module), "| Direction: ", unique(tmp$direction))
  p[[i]] <- ggplot(tmp, aes(x = reorder(ID, -p.adjust), 
                            y = -log10(p.adjust), 
                            fill = -log10(p.adjust))) +
    geom_bar(stat = "identity") +
    coord_flip() +
    xlab('') + ylab('−log10(adjusted p−value)') + ggtitle(title) +
    theme_bw() +
    theme_Publication(base_size = 12) 
}
ggsave(plot = wrap_plots(p, ncol = 1), 
       filename = file.path(patient_output_dir, 'ora_plots.pdf'), 
       width = 10, height = 10, device = 'pdf')
