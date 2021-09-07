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
                  fusion_method = fusion_method, 
                  snv_caller = snv_caller)

# run rna-seq analysis (output of this is required by several downstream scripts)
source(file.path(code_dir, "rnaseq_analysis", "run_rnaseq_analysis.R"))

# read sample information
sample_info <- read.delim(file.path(patient_dir, "clinical", "patient_report.txt"))
patient <- sample_info$subjectID

# update gsea enrichment output with each patient
if(snv_caller == "lancet"){
  gsea_enrichment <- file.path(code_dir, "rnaseq_analysis", 'gsea_enrichment.R')
  cmd <- paste('Rscript', gsea_enrichment, 
               '--patient', patient)
  print(cmd)
  system(cmd)
}

# generate genes and pathway enrichment output tables
if(snv_caller == "lancet"){
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

# key findings table
source(file.path(code_dir, "p1_modules", 'p1_key_clinical_findings.R'))

# disease specific information
source(file.path(code_dir, "p1_modules", 'p1_disease_specific_information.R'))  

# filter germline data
source(file.path(code_dir, "p1_modules", 'p1_filter_germline_vars.R'))              

# genomic summary table
source(file.path(code_dir, "p1_modules", 'p1_genomic_summary.R'))   

## page 2 modules 
# this is also part of the rnaseq_analysis module
# barplot of top 20 up/down genes
source(file.path(code_dir, "rnaseq_analysis", 'p2_diffexpr_genes_barplot.R'))

# barplot of top 10 up/down pathways
source(file.path(code_dir, "rnaseq_analysis", 'p2_diffreg_pathways_barplot.R'))  

## page 3 (tmb analysis)
# plot tumor mutational signatures
source(file.path(code_dir, "tmb_analysis", 'p3_tumor_signature_plot.R'))

# plot tumor mutational burden
source(file.path(code_dir, "tmb_analysis", 'p3_tmb_profile.R')) 

## prepare pbta + pnoc008 and tcga gbm + pnoc008 data for downstream functions 
if(snv_caller == "lancet"){
  source(file.path(code_dir, "transcriptomically_similar_analysis", 'tcga_format.R'))
  source(file.path(code_dir, "transcriptomically_similar_analysis", 'pbta_format.R'))
}

## page 4
# pediatric immune profiling using xcell (pbta)
source(file.path(code_dir, "immune_analysis", 'p4_immune_profile_pediatric.R'))

# adult immune profiling using xcell (tcga gbm)
source(file.path(code_dir, "immune_analysis", 'p4_immune_profile_adult.R'))

# tumor inflammation signature
source(file.path(code_dir, "immune_analysis", 'p4_tis_profile.R'))

## page 5 (pediatric analysis: pbta tumors)
# ssgsea
source(file.path(code_dir, "transcriptomically_similar_analysis", 'p5_ssgsea.R'))

# mutational analysis
source(file.path(code_dir, "transcriptomically_similar_analysis", 'p5_mutational_analysis_pediatric.R'))

# pathway analysis
source(file.path(code_dir, "transcriptomically_similar_analysis", 'p5_pathway_analysis_pediatric.R'))

# pediatric immune profiling using xcell (pbta transcriptomically correlated)
source(file.path(code_dir, "immune_analysis", 'p5_immune_profile_topcor_pediatric.R'))

# pediatric dimension reduction clustering
source(file.path(code_dir, "transcriptomically_similar_analysis", 'p5_dim_reduction_plot_pediatric.R'))

# pediatric km plot
source(file.path(code_dir, "survival_analysis", 'p5_kaplan_meier_pediatric.R'))

# pediatric transcriptomically similar patient table
source(file.path(code_dir, "transcriptomically_similar_analysis", 'p5_transcriptomically_similar_pediatric.R'))

## page 6 (adult analysis : TCGA tumors) 
# pathway analysis (top 20 transcriptomically similar patients)
source(file.path(code_dir, "transcriptomically_similar_analysis", 'p6_pathway_analysis_adult.R'))

# mutational analysis
source(file.path(code_dir, "transcriptomically_similar_analysis", 'p6_mutational_analysis_adult.R'))

# adult dimension reduction clustering
source(file.path(code_dir, "transcriptomically_similar_analysis", 'p6_dim_reduction_plot_adult.R'))

# adult km plot
source(file.path(code_dir, "survival_analysis", 'p6_kaplan_meier_adult.R'))

# adult transcriptomically similar patient table
source(file.path(code_dir, "transcriptomically_similar_analysis", 'p6_transcriptomically_similar_adult.R'))

## page 7
# genomic landscape plots
source(file.path(code_dir, "genomic_landscape_plots", "p7_circos_plot.R"))

# cnv plot (replaced with cnvkit's diagram.pdf)
# source(file.path(code_dir, "genomic_landscape_plots", "p7_cnv_plot.R"))

## page 8 (we don't need this currently)
# fname_cgs <- file.path(patient_dir, "output", "complexheatmap_cgs.png")
# fname_phgg <- file.path(patient_dir, "output", "complexheatmap_phgg.png")
# if(!file.exists(fname_cgs) | !file.exists(fname_phgg)){
#   source(file.path(patient_level_analyses, 'p8_cnv_exp_heatmap.R'))
# }

## page 9
# fname <- file.path(patient_dir, "output", "complexheatmap_oncogrid.pdf")
source(file.path(code_dir, "oncogrid_analysis", "p9_oncogrid.R"))

## page 10
# targeted findings
source(file.path(code_dir, "oncokb_analysis", 'p10_run_oncokb.R'))

# transcriptome based drug recommendations
source(file.path(code_dir, "drug_recommendations", "p10_transcriptome_drug_rec.R"))

# network plot - changes with snv caller (page 7)
# this is dependent on the output of transcriptome based drug recommendations
# so needs to be called after drug recommendations
source(file.path(code_dir, "genomic_landscape_plots", "p7_network_plot.R"))

# run CEMiTool to annotate hub genes in drug recommendations output
cemitools_script <- file.path(code_dir, "drug_recommendations", "p10_run_cemitools.R")
cmd <- paste('Rscript', cemitools_script, 
             '--patient', patient)
print(cmd)
system(cmd)

# dge density plots
source(file.path(code_dir, "drug_recommendations", "p10_drug_dge_density_plots.R"))

# drug pathways
source(file.path(code_dir, "drug_recommendations", 'p10_drug_pathways.R'))

# convert table outputs to visualizations
if(snv_caller == "lancet"){
  source(file.path(code_dir, "tables_to_plots", 'run_tables_to_plots.R'))
}

# run drug synergy module
if(snv_caller == "lancet"){
  source(file.path(code_dir, "drug_synergy", 'run_synergy.R'))
}

