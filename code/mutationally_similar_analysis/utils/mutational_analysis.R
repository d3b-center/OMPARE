# script to highlight relevant alterations top 20 transcriptomically similar patients

mutational_analysis <- function(nn_table, sample_id_interest, ref_cancer_dir, all_findings_output, key_clinical_findings_output, clinical, comparison){
  
  # tpm matrix of top 20 nearest neighbors
  top20 <- c(nn_table$nearest_neighbor, sample_id_interest)
  
  # subset clinical
  clinical <- list.files(path = ref_cancer_dir, pattern = "*histologies.tsv", full.names = T)
  clinical <- read.delim(clinical)
  clinical <- clinical %>%
    filter(sample_id %in% top20) %>%
    dplyr::select(Kids_First_Biospecimen_ID, sample_id, cohort, cohort_participant_id)
  
  # comparator data
  # mutations
  ref_maf <- list.files(path = ref_cancer_dir, pattern = "*mutation_filtered.rds", full.names = T)
  ref_maf <- readRDS(ref_maf)
  ref_maf <- ref_maf %>%
    filter(cohort_participant_id %in% clinical$cohort_participant_id) %>%
    mutate(Gene = Hugo_Symbol,
           Alteration_Datatype = "Mutation",
           Alteration_Type = Variant_Classification,
           Alteration = HGVSp_Short) %>%
    dplyr::select(Gene, Alteration_Datatype, Alteration_Type, Alteration, Kids_First_Biospecimen_ID, sample_id, cohort, cohort_participant_id) 

  # cnv 
  ref_cnv <- list.files(path = ref_cancer_dir, pattern = "*cnv_filtered.rds", full.names = T)
  ref_cnv <- readRDS(ref_cnv)
  ref_cnv <- ref_cnv %>%
    filter(cohort_participant_id %in% clinical$cohort_participant_id) %>%
    mutate(Gene = hgnc_symbol,
           Alteration_Datatype = "CNV",
           Alteration_Type = stringr::str_to_title(status),
           Alteration = paste0('Copy Number Value:', copy_number)) %>%
    dplyr::select(Gene, Alteration_Datatype, Alteration_Type, Alteration, Kids_First_Biospecimen_ID, sample_id, cohort, cohort_participant_id)
  
  # combine
  ref_variants <- rbind(ref_maf, ref_cnv)

  # add patient mutations
  patient_maf <- all_findings_output %>%
    ungroup(Type, Details) %>%
    filter(Type %in% c("Mutation", "VUS")) %>%
    mutate(Alteration_Datatype = "Mutation",
           Alteration_Type = gsub("Variant: | .*", "", Details),
           Alteration = gsub(".* ", "", Details)) %>%
    dplyr::rename("Gene" = "Aberration") %>%
    dplyr::select(colnames(ref_variants))
  
  # add patient copy number 
  patient_cnv <- all_findings_output %>%
    ungroup(Type, Details) %>%
    filter(Type %in% c("Amplification", "Loss", "Gain", "Complete Loss")) %>%
    mutate(Alteration_Datatype = "CNV",
           Alteration = gsub(" [|] .*", "", Details),
           Alteration_Type = gsub(".* ", "", Type)) %>%
    dplyr::rename("Gene" = "Aberration") %>%
    dplyr::select(colnames(ref_variants))
  
  # combine both
  patient_variants <- rbind(patient_maf, patient_cnv)
  
  # merge (all alterations from reference and from patient of interest)
  total_alterations <- unique(rbind(ref_variants, patient_variants))

  # alterations in genomically similar patients
  total_alt_table1 <- total_alterations %>%
    inner_join(total_alterations %>%
                 dplyr::select(Gene, Kids_First_Biospecimen_ID) %>%
                 unique() %>%
                 group_by(Gene) %>% 
                 dplyr::summarise(SampleCount = n()), by = c("Gene"))
  
  # at least 5/20 genomically similar patients
  total_alt_table1 <- total_alt_table1 %>%
    filter(SampleCount >= 5)
  
  # overlap with key clinical findings
  key_genes <- unique(key_clinical_findings_output$Aberration)
  total_alt_table2 <- total_alterations %>%
    filter(Gene %in% key_genes)
  
  # shared genes that are present in patient of interest + at least 1 more sample
  total_alt_table2 <- total_alt_table2 %>%
    inner_join(total_alt_table2 %>%
                 dplyr::select(Gene, Kids_First_Biospecimen_ID) %>% 
                 unique() %>%
                 group_by(Gene) %>% 
                 dplyr::summarise(SampleCount = n()), by = c("Gene")) %>%
    filter(SampleCount != 1)
  
  alt_tables <- list(recurrent_alterations = total_alt_table1, 
                     shared_genes = total_alt_table2)
  return(alt_tables)
}

