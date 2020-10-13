# Author: Komal S. Rathi
# Function: Up/Down genes + pathways for each PNOC008 sample and compare to GTEx Brain

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions for RNA-seq diffexpr & pathway analysis
source(file.path(patient_level_analyses_utils, "rnaseq_analysis_accessory.R"))

# format input expression data
# should be data-frame with 3 columns: Gene Sample TPM
expData.m <- expData %>%
  dplyr::rename('Gene' = 'gene_symbol',
                'TPM' = sampleInfo$subjectID) %>%
  mutate(Sample = sampleInfo$subjectID) %>%
  dplyr::select(Gene, Sample, TPM)
comparison = paste0(sampleInfo$subjectID, '_vs_GTExBrain')

# run single sample level analysis
rnaseq_analysis_output <- run_rnaseq_analysis(exp.data = expData.m,    # expression data (long format)
                                              refData = gtex_brain_tpm,      # gtex brain samples
                                              thresh = 2,              # foldchange cutoff 
                                              gene_set = gene_set, # reactome
                                              comparison = comparison, # comparison name
                                              single_sample = TRUE,    # single sample analysis
                                              sample_name = sampleInfo$subjectID) # pass sample name only for single sample analysis   

# save output
fname <- file.path(topDir, "output", "rnaseq_analysis_output.rds")
saveRDS(rnaseq_analysis_output, file = fname)