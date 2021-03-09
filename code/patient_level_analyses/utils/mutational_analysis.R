# script to highlight relevant alterations top 20 transcriptomically similar patients

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# reference directories
pbta_dir <- file.path(ref_dir, 'pbta')
tcga_dir <- file.path(ref_dir, 'tcga')
pnoc008_dir <- file.path(ref_dir, 'pnoc008')

mutational_analysis <- function(top_cor, key_clinical_findings_output, clinical, comparison){
  
  # matrix of top 20 correlated samples
  top20 <- colnames(top_cor)
  
  # subset clinical
  clinical <- clinical %>%
    filter(subject_id %in% top20) %>%
    dplyr::select(kids_first_biospecimen_id, sample_id, cohort_participant_id, subject_id)
  
  if(comparison == "pediatric"){
    # mutations, cnv, fusions
    tumor_mutations <- readRDS(file.path(pbta_dir, 'pbta-snv-consensus-mutation-filtered.rds'))
    tumor_cnv <- readRDS(file.path(pbta_dir, 'pbta-cnv-controlfreec-filtered.rds'))
    tumor_fusions <- readRDS(file.path(pbta_dir, 'pbta-fusion-putative-oncogenic-filtered.rds'))
  } else {
    # mutations and cnv
    tumor_mutations <- readRDS(file.path(tcga_dir, 'tcga_gbm_mutation_filtered.rds'))
    tumor_cnv <- readRDS(file.path(tcga_dir, 'tcga_gbm_cnv_filtered.rds'))
    tumor_fusions <- data.frame()
    clinical <- clinical %>%
      mutate(sample_id = ifelse(kids_first_biospecimen_id %in% c("", NA), subject_id, sample_id))
  }
  tumor_variants <- rbind(tumor_mutations, tumor_cnv, tumor_fusions)
  
  if(comparison == "adult"){
    tumor_variants$Kids_First_Biospecimen_ID <- ''
  }
  
  # combine with clinical (top 20)
  tumor_variants <- clinical %>%
    inner_join(tumor_variants, by = c("sample_id" = "SampleID")) %>%
    mutate(kids_first_biospecimen_id = Kids_First_Biospecimen_ID) %>%
    dplyr::select(-c(Kids_First_Biospecimen_ID))
  
  # merge pnoc008 mutations, copy number and fusions
  pnoc_mutations <- readRDS(file.path(pnoc008_dir, 'pnoc008_consensus_mutation_filtered.rds'))
  pnoc_cnv <- readRDS(file.path(pnoc008_dir, 'pnoc008_cnv_filtered.rds'))
  pnoc_fusions <- readRDS(file.path(pnoc008_dir, 'pnoc008_fusions_filtered.rds'))
  pnoc_variants <- rbind(pnoc_mutations, pnoc_cnv, pnoc_fusions)
  
  # combine with clinical (top 20)
  pnoc_variants <- clinical %>%
    inner_join(pnoc_variants, by = c("subject_id" = "SampleID")) %>%
    dplyr::select(-c(Kids_First_Biospecimen_ID))
  
  # merge
  total_alterations <- rbind(tumor_variants, pnoc_variants)
  
  # alterations in genomically similar patients
  total_alt_table1 <- total_alterations %>%
    inner_join(total_alterations %>%
                 dplyr::select(Gene, subject_id) %>%
                 unique() %>%
                 group_by(Gene) %>% 
                 summarise(SampleCount = n()), by = c("Gene"))
  
  # at least 5/20 genomically similar patients
  total_alt_table1 <- total_alt_table1 %>%
    filter(SampleCount >= 5)
  
  # overlap with key clinical findings
  key.clinical <- key_clinical_findings_output
  key.genes <- unique(key.clinical$Aberration)
  total_alt_table2 <- total_alterations %>%
    filter(Gene %in% key.genes)
  
  # shared genes that are present in patient of interest + at least 1 more sample
  total_alt_table2 <- total_alt_table2 %>%
    inner_join(total_alt_table2 %>%
                 dplyr::select(Gene, subject_id) %>% 
                 unique() %>%
                 group_by(Gene) %>% 
                 summarise(SampleCount = n()), by = c("Gene")) %>%
    filter(SampleCount != 1)
  
  alt_tables <- list(recurrent_alterations = total_alt_table1, 
                     shared_genes = total_alt_table2)
  return(alt_tables)
}
