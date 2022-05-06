# Author: Komal S. Rathi
# Function: driver script to call all functions and generate output
# NOTE: we will save the output of all modules which are independent of the somatic mutation caller

# set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
code_dir <- file.path(root_dir, "code")
utils_dir <- file.path(code_dir, "utils")

# load libraries
source(file.path(utils_dir, "load_libraries.R"))
`%>%` <- dplyr::`%>%`

# load reference data
source(file.path(utils_dir, 'load_reference.R')) 

# read patient data
source(file.path(utils_dir, "read_patient_data.R"))
read_patient_data(patient_dir = patient_dir, 
                  snv_caller = snv_caller)

# run rna-seq analysis (output of this is required by several downstream scripts)
output_dir <- file.path(patient_dir, "output", "rnaseq_analysis")
fname <- file.path(output_dir, "rnaseq_analysis_output.rds")
if(!file.exists(fname)){
  source(file.path(code_dir, "rnaseq_analysis", "run_rnaseq_analysis.R"))
} else {
  rnaseq_analysis_output <- readRDS(fname)
}

# read sample information
sample_info <- read.delim(file.path(patient_dir, "clinical", "patient_report.txt"))
patient <- sample_info$subjectID

# update gsea enrichment output with each patient
output_dir <- file.path(patient_dir, "output", "rnaseq_analysis")
fname <- file.path(output_dir, paste0(patient, "_summary_DE_Genes_Up.txt"))
if(!file.exists(fname)){
  # run gsea enrichment of patient vs other data sets
  gsea_enrichment <- file.path(code_dir, "rnaseq_analysis", 'gsea_enrichment.R')
  cmd <- paste('Rscript', gsea_enrichment, 
               '--patient', patient)
  print(cmd)
  system(cmd)
  
  # generate genes and pathway enrichment output tables
  enrichment_output <- file.path(code_dir, "rnaseq_analysis", 'enrichment_output.R')
  cmd <- paste('Rscript', enrichment_output, 
               '--patient', patient, 
               '--output', paste0(patient, '_summary'), 
               '--type', 'text')
  print(cmd)
  system(cmd)
}

## page 1 modules
# patient/sample information to display
source(file.path(code_dir, "p1_modules", 'p1_patient_sample_info.R')) 

# all findings table  
source(file.path(code_dir, "p1_modules", 'p1_all_findings.R'))

# key findings table (this is a subset of all_findings so will automatically get updated with above)
source(file.path(code_dir, "p1_modules", 'p1_key_clinical_findings.R'))

# disease specific information
source(file.path(code_dir, "p1_modules", 'p1_disease_specific_information.R'))  

# filter germline data
output_dir <- file.path(patient_dir, "output")
fname <- file.path(output_dir, "filtered_germline_vars.rds")
if(!file.exists(fname)){
  source(file.path(code_dir, "p1_modules", 'p1_filter_germline_vars.R'))              
} else {
  filtered_germ_vars <- readRDS(fname)
}

# genomic summary table
source(file.path(code_dir, "p1_modules", 'p1_genomic_summary.R'))   

## page 2 modules 
# this is also part of the rnaseq_analysis module
# barplot of top 20 up/down genes
output_dir <- file.path(patient_dir, "output", "rnaseq_analysis")
fname <- file.path(output_dir, "diffexpr_genes_barplot_output.pdf")
if(!file.exists(fname)){
  source(file.path(code_dir, "rnaseq_analysis", 'p2_diffexpr_genes_barplot.R'))
}

# barplot of top 10 up/down pathways
output_dir <- file.path(patient_dir, "output", "rnaseq_analysis")
fname <- file.path(output_dir, "diffreg_pathways_barplot_output.pdf")
if(!file.exists(fname)){
  source(file.path(code_dir, "rnaseq_analysis", 'p2_diffreg_pathways_barplot.R'))  
}

## page 3 (tmb analysis)
# plot tumor mutational signatures
output_dir <- file.path(patient_dir, "output", "tmb_analysis")
fname <- file.path(output_dir, "tumor_signature_output.rds")
if(!file.exists(fname)){
  source(file.path(code_dir, "tmb_analysis", 'p3_tumor_signature_plot.R'))
} else {
  tumor_signature_output <- readRDS(fname)
}

# plot tumor mutational burden
output_dir <- file.path(patient_dir, "output", "tmb_analysis")
fname <- file.path(output_dir, "tmb_profile_output.rds")
if(!file.exists(fname)){
  source(file.path(code_dir, "tmb_analysis", 'p3_tmb_profile.R')) 
} else {
  tmb_profile_output <- readRDS(fname)
}

## prepare tcga gbm + pnoc008 data for downstream functions 
if(snv_caller == "lancet"){
  source(file.path(code_dir, "transcriptomically_similar_analysis", 'tcga_format.R'))
} 

## prepare pbta + pnoc008 data for downstream functions
if(snv_caller == "lancet"){
  source(file.path(code_dir, "transcriptomically_similar_analysis", 'pbta_format.R'))
} 

## page 4
# pediatric immune profiling using xcell (pbta)
output_dir <- file.path(patient_dir, "output", "immune_analysis")
fname <- file.path(output_dir, 'immune_scores_pediatric.pdf')
if(!file.exists(fname)){
  source(file.path(code_dir, "immune_analysis", 'p4_immune_profile_pediatric.R'))
}

# adult immune profiling using xcell (tcga gbm)
output_dir <- file.path(patient_dir, "output", "immune_analysis")
fname <- file.path(output_dir, 'immune_scores_adult.pdf')
if(!file.exists(fname)){
  source(file.path(code_dir, "immune_analysis", 'p4_immune_profile_adult.R'))
}

# tumor inflammation signature
output_dir <- file.path(patient_dir, "output", "immune_analysis")
fname <- file.path(output_dir, 'tis_scores.pdf')
if(!file.exists(fname)){
  source(file.path(code_dir, "immune_analysis", 'p4_tis_profile.R'))
}

## page 5 (pediatric analysis: pbta tumors)
# ssgsea
output_dir <- file.path(patient_dir, "output", "transcriptomically_similar_analysis")
fname <- file.path(output_dir, 'ssgsea_scores_pediatric.pdf')
if(!file.exists(fname)){
  source(file.path(code_dir, "transcriptomically_similar_analysis", 'p5_ssgsea.R'))
}

# mutational analysis
output_dir <- file.path(patient_dir, "output", "transcriptomically_similar_analysis")
fname <- file.path(output_dir, "mutational_analysis_pediatric.rds")
if(!file.exists(fname)){
  source(file.path(code_dir, "transcriptomically_similar_analysis", 'p5_mutational_analysis_pediatric.R'))
}

# pathway analysis
output_dir <- file.path(patient_dir, "output", "transcriptomically_similar_analysis")
fname <- file.path(output_dir, "pathway_analysis_pediatric.rds")
if(!file.exists(fname)){
  source(file.path(code_dir, "transcriptomically_similar_analysis", 'p5_pathway_analysis_pediatric.R'))
}

# pediatric immune profiling using xcell (pbta transcriptomically correlated)
output_dir <- file.path(patient_dir, "output", "immune_analysis")
fname <- file.path(output_dir, 'immune_scores_topcor_pediatric.pdf')
if(!file.exists(fname)){
  source(file.path(code_dir, "immune_analysis", 'p5_immune_profile_topcor_pediatric.R'))
}

# pediatric dimension reduction clustering
output_dir <- file.path(patient_dir, "output", "transcriptomically_similar_analysis")
fname <- file.path(output_dir, "dim_reduction_plot_pediatric.rds")
if(!file.exists(fname)){
  source(file.path(code_dir, "transcriptomically_similar_analysis", 'p5_dim_reduction_plot_pediatric.R'))
} else {
  dim_reduction_plot_pediatric <- readRDS(fname)
}

# pediatric km plot
output_dir <- file.path(patient_dir, "output", "survival_analysis")
fname <- file.path(output_dir, "kaplan_meier_pediatric.pdf")
if(!file.exists(fname)){
  source(file.path(code_dir, "survival_analysis", 'p5_kaplan_meier_pediatric.R'))
} 

# pediatric transcriptomically similar patient table
output_dir <- file.path(patient_dir, "output", "transcriptomically_similar_analysis")
fname <- file.path(output_dir, "transciptomically_similar_pediatric.rds")
if(!file.exists(fname)){
  source(file.path(code_dir, "transcriptomically_similar_analysis", 'p5_transcriptomically_similar_pediatric.R'))
} else {
  transciptomically_similar_pediatric <- readRDS(fname)
}

## page 6 (adult analysis : TCGA tumors) 
# pathway analysis (top 20 transcriptomically similar patients)
output_dir <- file.path(patient_dir, "output", "transcriptomically_similar_analysis")
fname <- file.path(output_dir, "pathway_analysis_adult.rds")
if(!file.exists(fname)){
  source(file.path(code_dir, "transcriptomically_similar_analysis", 'p6_pathway_analysis_adult.R'))
}

# mutational analysis
output_dir <- file.path(patient_dir, "output", "transcriptomically_similar_analysis")
fname <- file.path(output_dir, "mutational_analysis_adult.rds")
if(!file.exists(fname)){
  source(file.path(code_dir, "transcriptomically_similar_analysis", 'p6_mutational_analysis_adult.R'))
}

# adult dimension reduction clustering
output_dir <- file.path(patient_dir, "output", "transcriptomically_similar_analysis")
fname <- file.path(output_dir, "dim_reduction_plot_adult.rds")
if(!file.exists(fname)){
  source(file.path(code_dir, "transcriptomically_similar_analysis", 'p6_dim_reduction_plot_adult.R'))
} else {
  dim_reduction_plot_adult <- readRDS(fname)
}

# adult km plot
output_dir <- file.path(patient_dir, "output", "survival_analysis")
fname <- file.path(output_dir, "kaplan_meier_adult.pdf")
if(!file.exists(fname)){
  source(file.path(code_dir, "survival_analysis", 'p6_kaplan_meier_adult.R'))
}

# adult transcriptomically similar patient table
output_dir <- file.path(patient_dir, "output", "transcriptomically_similar_analysis")
fname <- file.path(output_dir, "transciptomically_similar_adult.rds")
if(!file.exists(fname)){
  source(file.path(code_dir, "transcriptomically_similar_analysis", 'p6_transcriptomically_similar_adult.R'))
} else {
  transciptomically_similar_adult <- readRDS(fname)
}

## page 7
# genomic landscape plots
source(file.path(code_dir, "genomic_landscape_plots", "p7_circos_plot.R"))

# cnv plot (replaced with cnvkit's diagram.pdf)
# source(file.path(code_dir, "genomic_landscape_plots", "p7_cnv_plot.R"))

## page 8 (we don't need this currently)
# source(file.path(patient_level_analyses, 'p8_cnv_exp_heatmap.R'))

## page 9
output_dir <- file.path(patient_dir, "output", "oncogrid_analysis")
fname <- file.path(output_dir, "complexheatmap_oncogrid.pdf")
if(!file.exists(fname)){
  source(file.path(code_dir, "oncogrid_analysis", "p9_oncogrid.R"))
}

## page 10
# targeted findings
output_dir <- file.path(patient_dir, "output", "oncokb_analysis")
fname <- file.path(output_dir, paste0('oncokb_merged_', snv_caller, '_annotated_actgenes.txt'))
if(!file.exists(fname)){
  source(file.path(code_dir, "oncokb_analysis", 'p10_run_oncokb.R'))
} else {
  oncokb_output <- read.delim(fname, stringsAsFactors = F)
}

# update all findings with snv/indel hotspots and tier classifications
# this needs to be run after oncokb as it depends on its output
output_dir <- file.path(patient_dir, "output", "tier_classification")
fname <- file.path(output_dir, paste0("key_clinical_findings_output_", snv_caller, ".tsv"))
if(!file.exists(fname)){
  source(file.path(code_dir, 'tier_classification', 'run_tier_classification.R'))
}

# gsnca_module
# GSNCA analysis (required to annotate transcriptome_drug_rec.rds)
output_dir <- file.path(patient_dir, "output", "gsnca_analysis")
fname <- list.files(output_dir, "*.tsv")
if(length(fname) != 3){
  source(file.path(code_dir, "gsnca_analysis", "run_gsnca.R"))
}

# transcriptome based drug recommendations
output_dir <- file.path(patient_dir, "output", "drug_recommendations")
fname <- file.path(output_dir, "transcriptome_drug_rec.rds")
if(!file.exists(fname)){
  source(file.path(code_dir, "drug_recommendations", "p10_transcriptome_drug_rec.R"))
} else {
  transcriptome_drug_rec_output <- readRDS(fname)
}

# network plot - changes with snv caller (page 7)
# this is dependent on the output of transcriptome based drug recommendations
# so needs to be called after drug recommendations
source(file.path(code_dir, "genomic_landscape_plots", "p7_network_plot.R"))

# run CEMiTool to annotate hub genes in drug recommendations output
output_dir <- file.path(patient_dir, "output", "drug_recommendations")
fname <- file.path(output_dir, "ora_plots.pdf")
if(!file.exists(fname)){
  cemitools_script <- file.path(code_dir, "drug_recommendations", "p10_run_cemitools.R")
  cmd <- paste('Rscript', cemitools_script, 
               '--patient', patient)
  print(cmd)
  system(cmd)
}

# dge density plots
output_dir <- file.path(patient_dir, "output", "drug_recommendations", "drug_dge_density_plots")
fname <- file.path(output_dir, "top_drug_dge_density_plots.pdf")
if(!file.exists(fname)){
  source(file.path(code_dir, "drug_recommendations", "p10_drug_dge_density_plots.R"))
}

# drug pathways
output_dir <- file.path(patient_dir, "output", "drug_recommendations")
fname <- file.path(output_dir, "drug_pathways_barplot.pdf")
if(!file.exists(fname)){
  source(file.path(code_dir, "drug_recommendations", 'p10_drug_pathways.R'))
}

# convert table outputs to visualizations
if(snv_caller == "lancet"){
  source(file.path(code_dir, "tables_to_plots", 'run_tables_to_plots.R'))
}

# run drug synergy module
output_dir <- file.path(patient_dir, "output", "drug_synergy")
fname <- file.path(output_dir, "combined_qSig_synergy_score_top10.pdf")
if(!file.exists(fname)){
  source(file.path(code_dir, "drug_synergy", 'run_synergy.R'))
}

# Note: snv_caller == "lancet" condition is used for scripts that generate multiple outputs 
# which are not necessarily saved to files but are used for downstream analysis