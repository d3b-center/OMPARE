# script to highlight relevant alterations top 20 transcriptomically similar patients

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# reference directories
pbta_dir <- file.path(ref_dir, 'pbta')
tcga_dir <- file.path(ref_dir, 'tcga')
pnoc008_dir <- file.path(ref_dir, 'pnoc008')

mutational_analysis <- function(top_cor, key_clinical_findings_output, comparison){
  
  # matrix of top 20 correlated samples
  top20 <- colnames(top_cor)
  
  if(comparison == "pediatric"){
    # pbta corresponding sample ids
    pbta_clinical <- read.delim(file.path(pbta_dir, 'pbta-histologies.tsv'))
    tumor_clinical <- pbta_clinical %>%
      filter(Kids_First_Biospecimen_ID %in% top20)  %>%
      mutate(SampleID = sample_id) %>%
      dplyr::select(SampleID, Kids_First_Biospecimen_ID)
    
    # mutations, cnv, fusions
    tumor_mutations <- readRDS(file.path(pbta_dir, 'pbta-snv-consensus-mutation-filtered.rds'))
    tumor_cnv <- readRDS(file.path(pbta_dir, 'pbta-cnv-controlfreec-filtered.rds'))
    tumor_fusions <- readRDS(file.path(pbta_dir, 'pbta-fusion-putative-oncogenic-filtered.rds'))
    
  } else {
    # tcga clinical
    tcga_clinical <- readRDS(file.path(tcga_dir, 'tcga_gbm_clinical.rds'))
    tumor_clinical <- tcga_clinical %>%
      mutate(SampleID = sample_barcode,
             Kids_First_Biospecimen_ID = sample_barcode) %>%
      filter(Kids_First_Biospecimen_ID %in% top20) %>%
      dplyr::select(SampleID, Kids_First_Biospecimen_ID)
    
    # mutations and cnv
    tumor_mutations <- readRDS(file.path(tcga_dir, 'tcga_gbm_mutation_filtered.rds'))
    tumor_cnv <- readRDS(file.path(tcga_dir, 'tcga_gbm_cnv_filtered.rds'))
    tumor_fusions <- data.frame()
  }
  
  # pnoc008 corresponding sample ids
  pnoc008_clinical <- readRDS(file.path(pnoc008_dir, 'pnoc008_clinical.rds'))
  pnoc008_top20_clinical <- pnoc008_clinical %>%
    filter(subjectID %in% top20) %>%
    mutate(SampleID = subjectID, Kids_First_Biospecimen_ID = subjectID) %>%
    dplyr::select(SampleID, Kids_First_Biospecimen_ID)
  
  # combine other tumors with pnoc008
  combined_clinical <- rbind(tumor_clinical, pnoc008_top20_clinical)
  
  # merge other tumor + pnoc008 mutations, copy number and fusions
  # mutations
  pnoc_mutations <- readRDS(file.path(pnoc008_dir, 'pnoc008_consensus_mutation_filtered.rds'))
  total_mutations <- rbind(tumor_mutations, pnoc_mutations)
  
  # copy number
  pnoc_cnv <- readRDS(file.path(pnoc008_dir, 'pnoc008_cnv_filtered.rds'))
  total_cnv <- rbind(tumor_cnv, pnoc_cnv)
  
  # fusions
  pnoc_fusions <- readRDS(file.path(pnoc008_dir, 'pnoc008_fusions_filtered.rds'))
  total_fusions <- rbind(tumor_fusions, pnoc_fusions)
  
  # merge
  total_alterations <- rbind(total_mutations, total_cnv, total_fusions)
  
  # filter to top 20 genomically similar patients
  total_alterations <- total_alterations %>%
    filter(SampleID %in% combined_clinical$SampleID)
  
  # alterations in genomically similar patients
  total_alt_table1 <- total_alterations %>%
    inner_join(total_alterations %>%
                 dplyr::select(Gene, Kids_First_Biospecimen_ID) %>%
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
                 dplyr::select(Gene, SampleID) %>% 
                 unique() %>%
                 group_by(Gene) %>% 
                 summarise(SampleCount = n()), by = c("Gene")) %>%
    filter(SampleCount != 1)
  
  alt_tables <- list(recurrent_alterations = total_alt_table1, 
                     shared_genes = total_alt_table2)
  return(alt_tables)
}