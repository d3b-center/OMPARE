# Author: Komal S. Rathi
# Function: driver script to call all functions and generate output

# set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
code_dir <- file.path(root_dir, "code")
utils_dir <- file.path(code_dir, "utils")

# load libraries
source(file.path(utils_dir, "load_libraries.R"))
`%>%` <- dplyr::`%>%`

run_driver <- function(patient, patient_cancer_type, patient_dir){
  
  # directories
  normal_tissue_dir <- file.path(data_dir, "normal_data")
  pediatric_cancer_dir <- file.path(data_dir, "pediatric_data")
  pediatric_cancer_all_dir <- file.path(data_dir, "pediatric_data_all")
  adult_cancer_dir <- file.path(data_dir, "adult_data")
  
  # read patient data
  source(file.path(utils_dir, "read_patient_data.R"))
  read_patient_data(pediatric_cancer_dir = pediatric_cancer_dir,
                    patient_of_interest = patient,
                    mut_only = FALSE, 
                    rnaseq_only = FALSE)
  
  # run rna-seq analysis (output of this is required by several downstream scripts)
  output_dir <- file.path(patient_dir, "output", "rnaseq_analysis")
  fname <- file.path(output_dir, "rnaseq_analysis_output.rds")
  if(!file.exists(fname)){
    source(file.path(code_dir, "rnaseq_analysis", "run_rnaseq_analysis.R"))
  } else {
    rnaseq_analysis_output <- readRDS(fname)
  }
  
  # update gsea enrichment output with each patient
  source(file.path(code_dir, "rnaseq_analysis" ,"gsea_enrichment.R"))
  output_dir <- file.path(patient_dir, "output", "rnaseq_analysis")
  fname <- file.path(output_dir, "genes_down.txt")
  if(!file.exists(fname)){
    # run gsea enrichment of patient vs other data sets
    gsea_enrichment(normal_tissue = "Brain", 
                    adult_cancer = "GBM", 
                    pediatric_cancer = "HGAT", 
                    output_dir = output_dir,
                    tpm_data = tpm_data, 
                    count_data = count_data)
  }
  
  # all findings table  
  output_dir <- file.path(patient_dir, "output")
  fname <- file.path(output_dir, paste0("all_findings_output.rds"))
  if(!file.exists(fname)){
    source(file.path(code_dir, "p1_modules", 'p1_all_findings.R'))
  } else {
    all_findings_output <- readRDS(fname)
  }
  
  # key findings table 
  output_dir <- file.path(patient_dir, "output")
  fname <- file.path(output_dir, paste0("key_clinical_findings_output.rds"))
  if(!file.exists(fname)){
    source(file.path(code_dir, "p1_modules", 'p1_key_clinical_findings.R'))
  } else{
    key_clinical_findings_output <- readRDS(fname)
  }
  
  # disease specific information
  output_dir <- file.path(patient_dir, "output")
  fname <- file.path(output_dir, "disease_specific_information_output.rds")
  if(!file.exists(fname)){
    source(file.path(code_dir, "p1_modules", 'p1_disease_specific_information.R'))
  } else{
    disease_specific_information_output <- readRDS(fname)
  }
  
  # filter germline data
  output_dir <- file.path(patient_dir, "output")
  fname <- file.path(output_dir, "filtered_germline_vars.rds")
  if(!file.exists(fname)){
    source(file.path(code_dir, "p1_modules", 'p1_filter_germline_vars.R'))              
  } else {
    filtered_germ_vars <- readRDS(fname)
  }
  
  # genomic summary table
  output_dir <- file.path(patient_dir, "output")
  fname <- file.path(output_dir, "genomic_summary_output.rds")
  if(!file.exists(fname)){
    source(file.path(code_dir, "p1_modules", 'p1_genomic_summary.R'))
  } else {
    genomic_summary_output <- readRDS(fname)
  }
  
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
  
  # plot tumor mutational burden - needs mutect2 file
  output_dir <- file.path(patient_dir, "output", "tmb_analysis")
  fname <- file.path(output_dir, "tmb_profile_output.rds")
  if(!file.exists(fname)){
    source(file.path(code_dir, "tmb_analysis", 'p3_tmb_profile.R')) 
  } else {
    tmb_profile_output <- readRDS(fname)
  }
  
  # transcriptomically similar analysis 
  source(file.path(code_dir, "transcriptomically_similar_analysis", 'transcriptomically_similar_patients.R'))
  output_dir <- file.path(patient_dir, "output", "transcriptomically_similar_analysis")
  ## prepare pediatric tumors i.e. pbta + pnoc008 data for downstream functions
  fname <-  file.path(output_dir, paste0("pediatric_all_nn_table.rds"))
  if(!file.exists(fname)){
    tns_similar_analysis(ref_cancer_dir = pediatric_cancer_all_dir,
                         patient_dir = patient_dir,
                         sample_info = sample_info, 
                         tpm_data = tpm_data, 
                         prefix = "pediatric_all")
  } 
  
  ## prepare pediatric tumors i.e. pbta hgat + pnoc008 data for downstream functions
  fname <-  file.path(output_dir, paste0("pediatric_nn_table.rds"))
  if(!file.exists(fname)){
    tns_similar_analysis(ref_cancer_dir = pediatric_cancer_dir,
                         patient_dir = patient_dir,
                         sample_info = sample_info, 
                         tpm_data = tpm_data, 
                         prefix = "pediatric")
  } 
  
  # prepare adult tumors i.e. tcga + pnoc008 data for downstream functions
  fname <-  file.path(output_dir, paste0("adult_nn_table.rds"))
  if(!file.exists(fname)){
    tns_similar_analysis(ref_cancer_dir = adult_cancer_dir,
                         patient_dir = patient_dir,
                         sample_info = sample_info, 
                         tpm_data = tpm_data, 
                         prefix = "adult")
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
  # pediatric dimension reduction clustering
  output_dir <- file.path(patient_dir, "output", "transcriptomically_similar_analysis")
  fname <- file.path(output_dir, "dim_reduction_plot_pediatric.rds")
  if(!file.exists(fname)){
    source(file.path(code_dir, "transcriptomically_similar_analysis", 'p5_dim_reduction_plot_pediatric.R'))
  } else {
    dim_reduction_plot_pediatric <- readRDS(fname)
  }
  
  # pediatric transcriptomically similar patient table
  output_dir <- file.path(patient_dir, "output", "transcriptomically_similar_analysis")
  fname <- file.path(output_dir, "transciptomically_similar_pediatric.rds")
  if(!file.exists(fname)){
    source(file.path(code_dir, "transcriptomically_similar_analysis", 'p5_transcriptomically_similar_pediatric.R'))
  } else {
    transciptomically_similar_pediatric <- readRDS(fname)
  }
  
  # pediatric km plot
  output_dir <- file.path(patient_dir, "output", "survival_analysis")
  fname <- file.path(output_dir, "kaplan_meier_pediatric.pdf")
  if(!file.exists(fname)){
    source(file.path(code_dir, "survival_analysis", 'p5_kaplan_meier_pediatric.R'))
  } 
  
  # ssgsea
  output_dir <- file.path(patient_dir, "output", "transcriptomically_similar_analysis")
  fname <- file.path(output_dir, 'ssgsea_scores_pediatric.pdf')
  if(!file.exists(fname)){
    source(file.path(code_dir, "transcriptomically_similar_analysis", 'p5_ssgsea.R'))
  }
  
  # pediatric immune profiling using xcell (pbta transcriptomically correlated)
  output_dir <- file.path(patient_dir, "output", "immune_analysis")
  fname <- file.path(output_dir, 'immune_scores_topcor_pediatric.pdf')
  if(!file.exists(fname)){
    source(file.path(code_dir, "immune_analysis", 'p5_immune_profile_topcor_pediatric.R'))
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
  
  ## page 6 (adult analysis : TCGA tumors) 
  # adult dimension reduction clustering
  output_dir <- file.path(patient_dir, "output", "transcriptomically_similar_analysis")
  fname <- file.path(output_dir, "dim_reduction_plot_adult.rds")
  if(!file.exists(fname)){
    source(file.path(code_dir, "transcriptomically_similar_analysis", 'p6_dim_reduction_plot_adult.R'))
  } else {
    dim_reduction_plot_adult <- readRDS(fname)
  }
  
  # adult transcriptomically similar patient table
  output_dir <- file.path(patient_dir, "output", "transcriptomically_similar_analysis")
  fname <- file.path(output_dir, "transciptomically_similar_adult.rds")
  if(!file.exists(fname)){
    source(file.path(code_dir, "transcriptomically_similar_analysis", 'p6_transcriptomically_similar_adult.R'))
  } else {
    transciptomically_similar_adult <- readRDS(fname)
  }
  
  # adult km plot
  output_dir <- file.path(patient_dir, "output", "survival_analysis")
  fname <- file.path(output_dir, "kaplan_meier_adult.pdf")
  if(!file.exists(fname)){
    source(file.path(code_dir, "survival_analysis", 'p6_kaplan_meier_adult.R'))
  }
  
  # mutational analysis
  output_dir <- file.path(patient_dir, "output", "transcriptomically_similar_analysis")
  fname <- file.path(output_dir, "mutational_analysis_adult.rds")
  if(!file.exists(fname)){
    source(file.path(code_dir, "transcriptomically_similar_analysis", 'p6_mutational_analysis_adult.R'))
  }
  
  # pathway analysis (top 20 transcriptomically similar patients)
  output_dir <- file.path(patient_dir, "output", "transcriptomically_similar_analysis")
  fname <- file.path(output_dir, "pathway_analysis_adult.rds")
  if(!file.exists(fname)){
    source(file.path(code_dir, "transcriptomically_similar_analysis", 'p6_pathway_analysis_adult.R'))
  }
  
  ## page 7
  # genomic landscape plots
  output_dir <- file.path(patient_dir, "output", "genomic_landscape_plots")
  fname <- file.path(output_dir, "circos_plot.png")
  if(!file.exists(fname)){
    source(file.path(code_dir, "genomic_landscape_plots", "p7_circos_plot.R"))
  }
  
  ## page 8
  output_dir <- file.path(patient_dir, "output", "oncogrid_analysis")
  fname <- file.path(output_dir, "complexheatmap_oncogrid.pdf")
  if(!file.exists(fname)){
    source(file.path(code_dir, "oncogrid_analysis", "p8_oncogrid.R"))
  }
  
  ## page 9
  # targeted findings
  output_dir <- file.path(patient_dir, "output", "oncokb_analysis")
  fname <- file.path(output_dir, "oncokb_merged_consensus_annotated_actgenes.txt")
  if(!file.exists(fname)){
    source(file.path(code_dir, "oncokb_analysis", 'p9_run_oncokb.R'))
  } else {
    oncokb_output <- read.delim(fname, stringsAsFactors = F)
  }
  
  # update all findings with snv/indel hotspots and tier classifications
  # this needs to be run after oncokb as it depends on its output
  output_dir <- file.path(patient_dir, "output", "tier_classification")
  fname <- file.path(output_dir, "key_clinical_findings_output.tsv")
  if(!file.exists(fname)){
    source(file.path(code_dir, 'tier_classification', 'run_tier_classification.R'))
  }
  
  # clustering using coseq (only pediatric needed)
  output_dir <- file.path(patient_dir, "output", "coseq_detect", "pediatric")
  fname <- file.path(output_dir, "cancer_group_of_interest_nb_cluster_assigned.tsv")
  if(!file.exists(fname)){
    source(file.path(code_dir, "coseq_detect", 'run_coseq_detect_pediatric.R'))
  }
  
  # gsnca_module (only pediatric needed)
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
  
  # network plot (page 7)
  # this is dependent on the output of transcriptome based drug recommendations
  # so needs to be called after drug recommendations
  output_dir <- file.path(patient_dir, "output", "genomic_landscape_plots")
  fname <- file.path(output_dir, "network_plot_output.rds")
  if(!file.exists(fname)){
    source(file.path(code_dir, "genomic_landscape_plots", "p7_network_plot.R"))
  }
  
  # run CEMiTool to annotate hub genes in drug recommendations output
  output_dir <- file.path(patient_dir, "output", "drug_recommendations")
  fname <- file.path(output_dir, "ora_plots.pdf")
  if(!file.exists(fname)){
    source(file.path(code_dir, "drug_recommendations", "p10_run_cemitools.R"))
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
  
  # run drug synergy module
  output_dir <- file.path(patient_dir, "output", "drug_synergy")
  fname <- file.path(output_dir, "combined_qSig_synergy_score_top10.pdf")
  if(!file.exists(fname)){
    source(file.path(code_dir, "drug_synergy", 'run_synergy.R'))
  }
}
