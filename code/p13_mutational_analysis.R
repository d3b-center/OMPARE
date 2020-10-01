# Author: Komal S. Rathi
# Date: 06/29/2020
# Function: Script to highlight relevant mutational alterations
# top 20 genomically similar patients

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))
pbta.dir <- file.path(ref_dir, 'PBTA')
pnoc008.dir <- file.path(ref_dir, 'PNOC008')

recurrent.alterations <- function(topCor){
  
  # matrix of top 20 correlated samples
  top20 <- colnames(topCor)
  
  # PBTA corresponding sample ids
  pbta.clinData <- read.delim(file.path(pbta.dir, 'pbta-histologies.tsv'))
  pbta.clinData <- pbta.clinData %>%
    filter(Kids_First_Biospecimen_ID %in% top20)  %>%
    mutate(SampleID = sample_id) %>%
    dplyr::select(SampleID, Kids_First_Biospecimen_ID)
  
  # PNOC008 corresponding sample ids
  pnoc008.top20.clinData <- pnoc008.clinData %>%
    filter(subjectID %in% top20) %>%
    mutate(SampleID = subjectID, Kids_First_Biospecimen_ID = subjectID) %>%
    dplyr::select(SampleID, Kids_First_Biospecimen_ID)
  clinData <- rbind(pbta.clinData, pnoc008.top20.clinData)
  
  # merge PBTA + PNOC008 mutations, copy number and fusions
  
  # mutations
  pbta.mutations <- readRDS(file.path(pbta.dir, 'pbta-snv-consensus-mutation-filtered.rds'))
  pnoc.mutations <- readRDS(file.path(pnoc008.dir, 'PNOC008_consensus_mutData_filtered.rds'))
  total.mutations <- rbind(pbta.mutations, pnoc.mutations)
  
  # copy number
  pbta.cnv <- readRDS(file.path(pbta.dir, 'pbta-cnv-controlfreec-filtered.rds'))
  pnoc.cnv <- readRDS(file.path(pnoc008.dir, 'PNOC008_cnvData_filtered.rds'))
  total.cnv <- rbind(pbta.cnv, pnoc.cnv)
  
  # fusions
  pbta.fusions <- readRDS(file.path(pbta.dir, 'pbta-fusion-putative-oncogenic-filtered.rds'))
  pnoc.fusions <- readRDS(file.path(pnoc008.dir, 'PNOC008_fusData_filtered.rds'))
  total.fusions <- rbind(pbta.fusions, pnoc.fusions)
  
  # merge
  total.alterations <- rbind(total.mutations, total.cnv, total.fusions)
  
  # filter to top 20 genomically similar patients
  total.alterations <- total.alterations %>%
    filter(SampleID %in% clinData$SampleID)
  
  # alterations in genomically similar patients
  total.alt.table1 <- total.alterations %>%
    inner_join(total.alterations %>%
    dplyr::select(Gene, Kids_First_Biospecimen_ID) %>%
    unique() %>%
    group_by(Gene) %>% 
    summarise(SampleCount = n()), by = c("Gene"))
  
  # at least 5/20 genomically similar patients
  total.alt.table1 <- total.alt.table1 %>%
    filter(SampleCount >= 5)
  
  # overlap with key clinical findings
  key.clinical <- keyClinicalFindingsTable()
  key.genes <- unique(key.clinical$Aberration)
  total.alt.table2 <- total.alterations %>%
    filter(Gene %in% key.genes)
  
  # shared genes that are present in patient of interest + at least 1 more sample
  total.alt.table2 <- total.alt.table2 %>%
    inner_join(total.alt.table2 %>%
                 dplyr::select(Gene, SampleID) %>% 
                 unique() %>%
                 group_by(Gene) %>% 
                 summarise(SampleCount = n()), by = c("Gene")) %>%
    filter(SampleCount != 1)
  
  alt.tables <- list(total.alt.table1, total.alt.table2)
  return(alt.tables)
}

