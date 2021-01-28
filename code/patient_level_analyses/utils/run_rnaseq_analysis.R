# Author: Komal S. Rathi
# Function: Up/Down genes + pathways for each PNOC008 sample and compare to GTEx Brain

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions for RNA-seq diffexpr & pathway analysis
# source(file.path(patient_level_analyses_utils, "rnaseq_analysis_zscore.R"))
source(file.path(patient_level_analyses_utils, "rnaseq_analysis_edgeR.R"))

# format input expression data
# should be data-frame with 3 columns
# Gene Sample TPM
expData.m <- expData %>%
  dplyr::rename('Gene' = 'gene_symbol',
                'TPM' = sampleInfo$subjectID) %>%
  mutate(Sample = sampleInfo$subjectID) %>%
  dplyr::select(Gene, Sample, TPM)

# Gene Sample Counts
expData.counts.m <- expData.counts %>%
  dplyr::rename('Gene' = 'gene_symbol',
                'Counts' = sampleInfo$subjectID) %>%
  mutate(Sample = sampleInfo$subjectID) %>%
  dplyr::select(Gene, Sample, Counts)

# name of comparison
comparison = paste0(sampleInfo$subjectID, '_vs_GTExBrain')

# run single sample level z-score analysis
# rnaseq_analysis_output_zscore <- run_rnaseq_analysis(exp.data = expData.m, # tpm
#                                               refData = gtex_brain_tpm, # gtex tpm
#                                               thresh = 2, # foldchange cutoff 
#                                               gene_set = gene_set, # reactome
#                                               comparison = comparison, # comparison name
#                                               single_sample = TRUE,    # single sample analysis
#                                               sample_name = sampleInfo$subjectID) # pass sample name only for single sample analysis   

# run single sample level edgeR analysis
rnaseq_analysis_output_edger <- run_rnaseq_analysis_edger(exp.data.tpm = expData.m, # tpm
                                                          exp.data.counts = expData.counts.m, # counts
                                                          refData.counts = gtex_brain_counts, # gtex counts
                                                          gene_set = gene_set, # reactome
                                                          comparison = comparison, # comparison name
                                                          single_sample = TRUE,    # single sample analysis
                                                          sample_name = sampleInfo$subjectID) # pass sample name only for single sample analysis   
rnaseq_analysis_output <- rnaseq_analysis_output_edger

# save output
fname <- file.path(topDir, "output", "rnaseq_analysis_output.rds")
saveRDS(rnaseq_analysis_output, file = fname)
