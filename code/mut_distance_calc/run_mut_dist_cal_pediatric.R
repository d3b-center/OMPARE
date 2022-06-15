# define directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "mut_distance_calc")
output_dir <- file.path(patient_dir, "output", "mut_distance_calc", "pediatric")
dir.create(output_dir, showWarnings = F, recursive = T)

# sample_info is sample info on patient of interest
sample_id_interest <- sample_info %>%
  pull(sample_id) %>% 
  unique()

# other references file
genes_df <- file.path(data_dir, "oncokb", "OncoKB_DGD_comb_genes.txt")

# generate binary matrix
binary_mat_file <- file.path(output_dir, "binary_mutational_matrix.tsv")
make_binary_script <- file.path(module_dir, '01-make_binary_matrix.R')
cmd <- paste('Rscript', make_binary_script,
             '--ref_cancer_dir', pediatric_cancer_dir, 
             '--mat_file', binary_mat_file,
             '--genes', genes_df,
             '--sample_id_interest', sample_id_interest,
             '--analysis_type', 'pediatric')
system(cmd)

# calculate mutational distance
distance_matrix <- file.path(output_dir, "poi_vs_cohort_distance.tsv")
calc_dist_script <- file.path(module_dir, '02-calc_mutational_dist.R')
cmd <- paste('Rscript', calc_dist_script,
             '--sample_id_interest', sample_id_interest,
             '--mat_file', binary_mat_file,
             '--dist_file', distance_matrix)
system(cmd)

# determine the cluster assignment of the 20 closest cases in mutational space 
coseq_cluster_assignment <- file.path(patient_dir, "output", "coseq_detect", "pediatric", "cancer_group_of_interest_nb_cluster_assigned.tsv")
cluster_by_mut_dist_script <- file.path(module_dir, '03-cluster_by_mut_dist.R')
cmd <- paste('Rscript', cluster_by_mut_dist_script ,
             '--ref_cancer_dir', pediatric_cancer_dir, 
             '--sample_id_interest', sample_id_interest,
             '--coseq_cluster_assignment', coseq_cluster_assignment,
             '--mut_dist', distance_matrix,
             '--output_dir', output_dir)
system(cmd)

# calculate the v-test score for the cluster of interest
vtest_file <- file.path(output_dir, "cluster_gene_vscore.rds")
filtered_count_cg_coding_pval <- file.path(patient_dir, "output", "coseq_detect", "pediatric", "filtered_count_cg_coding_pval.rds")
vtest_score_calc_script <- file.path(module_dir, '04-vtest_scores_cal.R')
cmd <- paste('Rscript', vtest_score_calc_script ,
             '--ref_cancer_dir', pediatric_cancer_dir, 
             '--coseq_cluster_assignment', coseq_cluster_assignment,
             '--filtered_count_cg_coding_pval', filtered_count_cg_coding_pval,
             '--vtest_file', vtest_file)
system(cmd)

# run drug signature analysis using the Gene and V-test score as input
nb_cluster_assigned <- file.path(patient_dir, "output", "mut_distance_calc", "pediatric", "patient_of_interest_nb_cluster_assigned.tsv")
linc_connectivity_script <- file.path(module_dir, '05-linc_based_connectivity.R')
cmd <- paste('Rscript', linc_connectivity_script ,
             '--patient', sample_id_interest,
             '--nb_cluster_assigned', nb_cluster_assigned,
             '--vtest_output', vtest_file, 
             '--wtcs_fdr_cutoff', 0.05,
             '--trend_val', 'down',
             '--cor_score_cutoff', 0, 
             '--num_sets', 25,
             '--output_dir', output_dir)
system(cmd)
