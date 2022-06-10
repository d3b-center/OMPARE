# Author: Komal S. Rathi
# Function: driver script to call all functions and generate output for mutation only analysis

# set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
code_dir <- file.path(root_dir, "code")
utils_dir <- file.path(code_dir, "utils")

# load libraries
source(file.path(utils_dir, "load_libraries.R"))
`%>%` <- dplyr::`%>%`

run_driver <- function(patient, patient_cancer_type, snv_caller, patient_dir){
  
  # directories
  normal_tissue_dir <- file.path(data_dir, "normal_data")
  pediatric_cancer_dir <- file.path(data_dir, "pediatric_data")
  adult_cancer_dir <- file.path(data_dir, "adult_data")
  
  # read patient data
  source(file.path(utils_dir, "read_patient_data.R"))
  read_patient_data(pediatric_cancer_dir = pediatric_cancer_dir, 
                    patient_of_interest = patient,
                    mut_only = TRUE,
                    rnaseq_only = FALSE,
                    snv_caller = snv_caller)
  
  # assign NULL to rnaseq_analysis_output
  rnaseq_analysis_output <- NULL
  
  # all findings table  
  output_dir <- file.path(patient_dir, "output")
  fname <- file.path(output_dir, "all_findings_output.rds")
  if(!file.exists(fname)){
    source(file.path(code_dir, "p1_modules", 'p1_all_findings.R'))
  } else {
    all_findings_output <- readRDS(fname)
  }
  
  # key findings table 
  output_dir <- file.path(patient_dir, "output")
  fname <- file.path(output_dir, "key_clinical_findings_output.rds")
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
  
  # filter germline data - for germline file, download patient specific file from cavatica
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
  
  # clustering using coseq (pediatric)
  output_dir <- file.path(patient_dir, "output", "coseq_detect", "pediatric")
  fname <- file.path(output_dir, "cancer_group_of_interest_nb_cluster_assigned.tsv")
  if(!file.exists(fname)){
    source(file.path(code_dir, "coseq_detect", 'run_coseq_detect_pediatric.R'))
  }
  
  # mutational distance calculation (pediatric)
  # drug pathways
  output_dir <- file.path(patient_dir, "output", "mut_distance_calc", "pediatric")
  fname <- file.path(output_dir, "drug_pathways_barplot.pdf")
  if(!file.exists(fname)){
    source(file.path(code_dir, "mut_distance_calc", 'run_mut_dist_cal_pediatric.R'))
  }
  
  ## mutationally similar analysis
  
  # survival analysis will use the above info from coseq+mut distance calculation
  # km plot (pediatric)
  output_dir <- file.path(patient_dir, "output", "mutationally_similar_analysis")
  fname <- file.path(output_dir, "kaplan_meier_pediatric.pdf")
  if(!file.exists(fname)){
    source(file.path(code_dir, "mutationally_similar_analysis", 'kaplan_meier_pediatric.R'))
  }
  
  # mutational analysis (pediatric)
  output_dir <- file.path(patient_dir, "output", "mutationally_similar_analysis")
  fname <- file.path(output_dir, "mutational_analysis_pediatric.rds")
  if(!file.exists(fname)){
    source(file.path(code_dir, "mutationally_similar_analysis", 'mutational_analysis_pediatric.R'))
  }
  
  # clustering using coseq (adult tumors)
  output_dir <- file.path(patient_dir, "output", "coseq_detect", "adult")
  fname <- file.path(output_dir, "cancer_group_of_interest_nb_cluster_assigned.tsv")
  if(!file.exists(fname)){
    source(file.path(code_dir, "coseq_detect", 'run_coseq_detect_adult.R'))
  }
  
  # mutational distance calculation (adult tumors)
  # drug pathways
  output_dir <- file.path(patient_dir, "output", "mut_distance_calc", "adult")
  fname <- file.path(output_dir, "drug_pathways_barplot.pdf")
  if(!file.exists(fname)){
    source(file.path(code_dir, "mut_distance_calc", 'run_mut_dist_cal_adult.R'))
  }
  
  ## mutationally similar analysis
  
  # survival analysis will use the above info from coseq+mut distance calculation
  # km plot (adult)
  output_dir <- file.path(patient_dir, "output", "mutationally_similar_analysis")
  fname <- file.path(output_dir, "kaplan_meier_adult.pdf")
  if(!file.exists(fname)){
    source(file.path(code_dir, "mutationally_similar_analysis", 'kaplan_meier_adult.R'))
  }
  
  # mutational analysis (adult)
  output_dir <- file.path(patient_dir, "output", "mutationally_similar_analysis")
  fname <- file.path(output_dir, "mutational_analysis_adult.rds")
  if(!file.exists(fname)){
    source(file.path(code_dir, "mutationally_similar_analysis", 'mutational_analysis_adult.R'))
  }
  
  # oncogrid
  output_dir <- file.path(patient_dir, "output", "oncogrid_analysis")
  fname <- file.path(output_dir, "complexheatmap_oncogrid.pdf")
  if(!file.exists(fname)){
    source(file.path(code_dir, "oncogrid_analysis", "p8_oncogrid.R"))
  }
  
  # targeted findings (uses unfiltered cnv, fusions and snv)
  output_dir <- file.path(patient_dir, "output", "oncokb_analysis")
  fname <- file.path(output_dir, paste0('oncokb_merged_', snv_caller, '_annotated_actgenes.txt'))
  if(!file.exists(fname)){
    source(file.path(code_dir, "oncokb_analysis", 'p9_run_oncokb.R'))
  } else {
    oncokb_output <- read.delim(fname, stringsAsFactors = F)
  }
  
  # update all findings and key findings with snv/indel hotspots and tier classifications
  # this needs to be run after oncokb as it depends on its output
  output_dir <- file.path(patient_dir, "output", "tier_classification")
  fname <- file.path(output_dir, paste0("key_clinical_findings_output_", snv_caller, ".tsv"))
  if(!file.exists(fname)){
    source(file.path(code_dir, 'tier_classification', 'run_tier_classification.R'))
  }
  
  # run CEMiTool to annotate hub genes in drug recommendations output (pediatric)
  output_dir <- file.path(patient_dir, "output", "drug_recommendations", "pediatric")
  fname <- file.path(output_dir, "ora_plots.pdf")
  if(!file.exists(fname)){
    source(file.path(code_dir, "drug_recommendations", "run_cemitools_mut_only_pediatric.R"))
  }
  
  # run CEMiTool to annotate hub genes in drug recommendations output (adult)
  output_dir <- file.path(patient_dir, "output", "drug_recommendations", "adult")
  fname <- file.path(output_dir, "ora_plots.pdf")
  if(!file.exists(fname)){
    source(file.path(code_dir, "drug_recommendations", "run_cemitools_mut_only_adult.R"))
  }
}