# Author: Komal S. Rathi
# Function: driver script to call all functions and generate output
# NOTE: we will save the output of all modules which are independent of the somatic mutation caller

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# load reference data
source(file.path(utils_dir, 'load_reference.R')) 

# read patient data
source(file.path(patient_level_analyses_utils, "read_patient_data.R"))
readData(topDir = topDir, fusion_method = fusion_method, snv_caller = snv_caller)

# run rna-seq analysis (output of this is required by several downstream scripts)
fname <- file.path(topDir, "output", "rnaseq_analysis_output.rds")
if(file.exists(fname)){
  rnaseq_analysis_output <- readRDS(fname)
} else {
  source(file.path(patient_level_analyses_utils, "run_rnaseq_analysis.R"))
}

## page 1 modules
# patient/sample information to display
source(file.path(patient_level_analyses, 'p1_patient_sample_info.R'))    

# all findings table  
source(file.path(patient_level_analyses, 'p1_all_findings.R'))

# key findings table
source(file.path(patient_level_analyses, 'p1_key_clinical_findings.R'))

# disease specific information
source(file.path(patient_level_analyses, 'p1_disease_specific_information.R'))  

# filter germline data
fname <- file.path(topDir, "output", "filtered_germline_vars.rds")
if(file.exists(fname)){
  filtered_germ_vars <- readRDS(fname)
} else {
  source(file.path(patient_level_analyses, 'p1_filter_germline_vars.R'))              
}

# genomic summary table
source(file.path(patient_level_analyses, 'p1_genomic_summary.R'))   

## page 2 modules
# barplot of top 20 up/down genes
fname <- file.path(topDir, "output", "diffexpr_genes_barplot_output.rds")
if(file.exists(fname)){
  diffexpr_genes_barplot_output <- readRDS(fname)
} else {
  source(file.path(patient_level_analyses, 'p2_diffexpr_genes_barplot.R'))              
}

# barplot of top 10 up/down pathways
fname <- file.path(topDir, "output", "diffreg_pathways_barplot_output.rds")
if(file.exists(fname)){
  diffreg_pathways_barplot_output <- readRDS(fname)
} else {
  source(file.path(patient_level_analyses, 'p2_diffreg_pathways_barplot.R'))              
}

## page 3
# plot tumor mutational signatures
fname <- file.path(topDir, "output", "tumor_signature_output.rds")
if(file.exists(fname)){
  tumor_signature_output <- readRDS(fname)
} else {
  source(file.path(patient_level_analyses, 'p3_tumor_signature_plot.R'))
}

# plot tumor mutational burden
fname <- file.path(topDir, "output", "tmb_profile_output.rds")
if(file.exists(fname)){
  tmb_profile_output <- readRDS(fname)
} else {
  source(file.path(patient_level_analyses, 'p3_tmb_profile.R')) 
}

# prepare pbta and tcga gbm data for downstream functions 
source(file.path(patient_level_analyses, 'tcga_format.R'))
source(file.path(patient_level_analyses, 'pbta_format.R'))

## page 4
# pediatric immune profiling using xcell (pbta)
source(file.path(patient_level_analyses, 'p4_immune_profile_pediatric.R'))

# adult immune profiling using xcell (tcga gbm)
source(file.path(patient_level_analyses, 'p4_immune_profile_adult.R'))

# tumor inflammation signature
source(file.path(patient_level_analyses, 'p4_tis_profile.R'))

## page 5 (pediatric analysis: pbta tumors)
# ssgsea
source(file.path(patient_level_analyses, 'p5_ssgsea.R'))

# mutational analysis
fname <- file.path(topDir, "output", "mutational_analysis_pediatric.rds")
if(file.exists(fname)){
  mutational_analysis_pediatric <- readRDS(fname)
} else {
  source(file.path(patient_level_analyses, 'p5_mutational_analysis_pediatric.R'))
}

# pathway analysis
fname <- file.path(topDir, "output", "pathway_analysis_pediatric.rds")
if(file.exists(fname)){
  pathway_analysis_pediatric <- readRDS(fname)
} else {
  source(file.path(patient_level_analyses, 'p5_pathway_analysis_pediatric.R'))
}

# pediatric immune profiling using xcell (pbta transcriptomically correlated)
source(file.path(patient_level_analyses, 'p5_immune_profile_topcor_pediatric.R'))

# pediatric dimention reduction clustering
fname <- file.path(topDir, "output", "dim_reduction_plot_pediatric.rds")
if(file.exists(fname)){
  dim_reduction_plot_pediatric <- readRDS(fname)
} else {
  source(file.path(patient_level_analyses, 'p5_dim_reduction_plot_pediatric.R'))
}

# pediatric km plot
fname <- file.path(topDir, "output", "kaplan_meier_pediatric.rds")
if(file.exists(fname)){
  kaplan_meier_pediatric <- readRDS(fname)
} else {
  source(file.path(patient_level_analyses, 'p5_kaplan_meier_pediatric.R'))
}

# pediatric transcriptomically similar patient table
fname <- file.path(topDir, "output", "transciptomically_similar_pediatric.rds")
if(file.exists(fname)){
  transciptomically_similar_pediatric <- readRDS(fname)
} else {
  source(file.path(patient_level_analyses, 'p5_transcriptomically_similar_pediatric.R'))
}

## page 6 (adult analysis : TCGA tumors) 
# pathway analysis (top 20 transcriptomically similar patients)
fname <- file.path(topDir, "output", "pathway_analysis_adult.rds")
if(file.exists(fname)){
  pathway_analysis_adult <- readRDS(fname)
} else {
  source(file.path(patient_level_analyses, 'p6_pathway_analysis_adult.R'))
}

# mutational analysis
fname <- file.path(topDir, "output", "mutational_analysis_adult.rds")
if(file.exists(fname)){
  mutational_analysis_adult <- readRDS(fname)
} else {
  source(file.path(patient_level_analyses, 'p6_mutational_analysis_adult.R'))
}

# adult dimension reduction clustering
fname <- file.path(topDir, "output", "dim_reduction_plot_adult.rds")
if(file.exists(fname)){
  dim_reduction_plot_adult <- readRDS(fname)
} else {
  source(file.path(patient_level_analyses, 'p6_dim_reduction_plot_adult.R'))
}

# adult km plot
fname <- file.path(topDir, "output", "kaplan_meier_adult.rds")
if(file.exists(fname)){
  kaplan_meier_adult <- readRDS(fname)
} else {
  source(file.path(patient_level_analyses, 'p6_kaplan_meier_adult.R'))
}

# adult transcriptomically similar patient table
fname <- file.path(topDir, "output", "transciptomically_similar_adult.rds")
if(file.exists(fname)){
  transciptomically_similar_adult <- readRDS(fname)
} else {
  source(file.path(patient_level_analyses, 'p6_transcriptomically_similar_adult.R'))
}

## page 7
source(file.path(patient_level_analyses, 'p7_circos_plot.R'))

# cnv plot
fname <- file.path(topDir, "output", "cnv_plot.png")
if(!file.exists(fname)){
  source(file.path(patient_level_analyses, 'p7_cnv_plot.R'))
}

## page 8 
# fname_cgs <- file.path(topDir, "output", "complexheatmap_cgs.png")
# fname_phgg <- file.path(topDir, "output", "complexheatmap_phgg.png")
# if(!file.exists(fname_cgs) | !file.exists(fname_phgg)){
#   source(file.path(patient_level_analyses, 'p8_cnv_exp_heatmap.R'))
# }

## page 9
fname <- file.path(topDir, "output", "complexheatmap_oncogrid.pdf")
if(!file.exists(fname)){
  source(file.path(patient_level_analyses, 'p9_oncogrid.R'))
}

## page 10
# targeted findings
fname <- file.path(topDir, "output", paste0('oncokb_merged_', snv_caller, '_annotated_actgenes.txt'))
if(file.exists(fname)){
  oncokb_output <- read.delim(fname, stringsAsFactors = F)
} else {
  source(file.path(patient_level_analyses, 'p10_run_oncokb.R'))
  oncokb_output <- read.delim(fname, stringsAsFactors = F)
}

# transcriptome based drug recommendations
fname <- file.path(topDir, "output", "transcriptome_drug_rec.rds")
if(file.exists(fname)){
  transcriptome_drug_rec_output <- readRDS(fname)
} else {
  source(file.path(patient_level_analyses, 'p10_transcriptome_drug_rec.R'))
}

# network plot - changes with snv caller (page 7)
# this is dependent on the output of transcriptome based drug recommendations
source(file.path(patient_level_analyses, 'p7_network_plot.R'))

# run CEMiTool to annotate hub genes on above output
test_file <- file.path(topDir, "output", "ora_plots.pdf")
if(file.exists(test_file)){
  fname <- file.path(topDir, "output", "transcriptome_drug_rec.rds")
  transcriptome_drug_rec_output <- readRDS(fname)
} else {
  source(file.path(patient_level_analyses, 'run_cemitools.R'))
}

# dge density plots
fname <- file.path(topDir, 'output', 'drug_dge_density_plots', 'top_drug_dge_density_plots.pdf')
if(!file.exists(fname)){
  source(file.path(patient_level_analyses, 'p10_drug_dge_density_plots.R'))
}

# drug pathways
fname <- file.path(topDir, "output", "drug_pathways_barplot.pdf")
if(!file.exists(fname)){
  source(file.path(patient_level_analyses, 'p10_drug_pathways.R'))
}

# convert table outputs to visualizations
source(file.path(patient_level_analyses, 'run_tables_to_plots.R'))
