# Author: Komal S. Rathi
# Function: Up/Down genes + pathways for each PNOC008 sample and compare to GTEx Brain

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "rnaseq_analysis")
output_dir <- file.path(patient_dir, "output", "rnaseq_analysis")
dir.create(output_dir, recursive = T, showWarnings = F)

# source functions for RNA-seq diffexpr & pathway analysis
source(file.path(module_dir, "utils", "rnaseq_analysis_edgeR.R"))

# cancer genes
cancer_genes <- readRDS(file.path(root_dir, "data", "cancer_gene_list.rds"))

# format input expression data
expData.m <- expData %>% 
  dplyr::select(-c(gene_id)) %>%
  gather(sample, tpm, -c('gene_symbol'))

# Gene Sample Counts
expData.counts.m <- expData.counts %>% 
  dplyr::select(-c(gene_id)) %>%
  gather(sample, expected_count, -c('gene_symbol'))

# name of comparison
comparison = paste0(sampleInfo$subjectID, '_vs_GTExBrain')

# module directory
output_dir <- file.path(patient_dir, "output", "rnaseq_analysis")
dir.create(output_dir, recursive = T, showWarnings = F)

# output file
fname <- file.path(output_dir, "rnaseq_analysis_output.rds")

if(!file.exists(fname)){
  # run single sample level edgeR analysis
  rnaseq_analysis_output <- run_rnaseq_analysis_edger(exp.data.tpm = expData.m, # tpm
                                                      exp.data.counts = expData.counts.m, # counts
                                                      refData.counts = gtex_brain_counts, # gtex counts
                                                      gene_set = gene_set, # reactome
                                                      comparison = comparison, # comparison name
                                                      single_sample = TRUE,    # single sample analysis
                                                      sample_name = sampleInfo$subjectID, # pass sample name only for single sample analysis 
                                                      cancer_genes = cancer_genes)   
  
  # save output
  saveRDS(rnaseq_analysis_output, file = fname)
} else {
  print("output exists")
  rnaseq_analysis_output <- readRDS(fname)
}