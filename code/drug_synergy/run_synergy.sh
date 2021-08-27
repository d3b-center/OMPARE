#!/bin/bash

# Author: Run Jin

# This analysis calculates the drug synergy score with the workflow laid-out in this ticket
# (https://github.com/d3b-center/bixu-tracker/issues/1202)

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
analysis_dir="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$analysis_dir" || exit

data_dir="../../data"
results_dir="results"
results_dir_synergy="${results_dir}/synergy_score"
results_dir_map="${results_dir}/drug_gene_map"
ref_dir="../../references"

# Files for generating subnetworks
interaction_file="${data_dir}/interactions.tsv"
enrichment_nes_file="${data_dir}/enrichment_nes.tsv"
cluster_file="${data_dir}/CC_based_heatmap_Distance_euclidean_finalLinkage_average_clusterAlg_KM_expct_counts_VST_cluster_and_annotation.tsv"

# qSig files from CEMiTools output
gtex_qSig="${data_dir}/GTExBrain_qSig_output.txt"
pbta_qSig="${data_dir}/PBTA_ALL_qSig_output.txt"
pbta_hgg_qSig="${data_dir}/PBTA_HGG_qSig_output.txt"

# subnetwork file and drug mapped subnetwork qSig file 
subnetwork="${results_dir}/subnetwork_genes.tsv"
subnetwork_mapped="${results_dir_map}/subnetwork_gene_drug_map.tsv"

gtex_subnet_qSig_mapped="${results_dir_map}/gtex_qSig_subnetwork_drug_gene_map.tsv"
pbta_subnet_qSig_mapped="${results_dir_map}/pbta_qSig_subnetwork_drug_gene_map.tsv"
pbta_hgg_subnet_qSig_mapped="${results_dir_map}/pbta_hgg_qSig_subnetwork_drug_gene_map.tsv"

# References database database path
chembldb_path="${ref_dir}/chembl_29_sqlite/chembl_29.db"

# Obtain drugs that are associated with all the genes in the subnetwork
Rscript --vanilla ${analysis_dir}/01-subnetwork_qSig_gene_drug_map.R \
--interaction $interaction_file \
--enrichment_nes $enrichment_nes_file \
--cluster $cluster_file \
--chemblDb_path $chembldb_path \
--sample_interest "P-38" \
--gtex_qSig $gtex_qSig \
--pbta_qSig $pbta_qSig \
--pbta_hgg_qSig $pbta_hgg_qSig \
--subnetwork $subnetwork \
--subnetwork_mapped $subnetwork_mapped \
--gtex_mapped $gtex_subnet_qSig_mapped \
--pbta_mapped $pbta_subnet_qSig_mapped \
--pbta_hgg_mapped $pbta_hgg_subnet_qSig_mapped

# Obtain drugs that are both in qSig and subnetwork
Rscript --vanilla ${analysis_dir}/02-drug_synergy_score_calc.R \
--gtex_mapped $gtex_subnet_qSig_mapped \
--pbta_mapped $pbta_subnet_qSig_mapped \
--pbta_hgg_mapped $pbta_hgg_subnet_qSig_mapped \
--subnetwork $subnetwork \
--output_path $results_dir_synergy




