# function: high confidence findings
suppressPackageStartupMessages({
  library(tidyverse)
})

# gene sets (c2 reactome)
gene_set <- getGmt(file.path(data_dir, 'msigdb', 'c2.cp.reactome.v6.0.symbols.gmt'), collectionType = BroadCollection(), geneIdType = SymbolIdentifier())
gene_set <- geneIds(gene_set)
gene_set_ts <- stack(gene_set)

key_clinical_findings <- function(all_findings_output, rnaseq_analysis_output) {

  # use all findings as base
  myTable <- all_findings_output
  
  # annotate with expression data
  if(!is.null(rnaseq_analysis_output)){
    
    # gene expression as tpm
    rnaEvidence <- rnaseq_analysis_output$tpm %>%
      dplyr::rename(tpm = 1) %>%
      rownames_to_column("Gene")
    
    # map genes to significant pathways
    sigGeneSets <- rnaseq_analysis_output$pathways %>%
      dplyr::select(pathway, direction) %>%
      inner_join(gene_set_ts, by = c("pathway" = "ind"))  %>%
      dplyr::mutate(pathway = paste0(pathway, "(", direction, ")")) %>%
      dplyr::select(pathway, values) %>%
      group_by(values) %>%
      dplyr::summarise(Pathway = toString(pathway))
    
    # filter to deletions that have TPM < 10
    myTableDel <- myTable %>%
      filter(Type %in% c("Loss", "Complete Loss")) %>%
      left_join(rnaEvidence, by = c("Aberration" = "Gene")) %>%
      left_join(sigGeneSets, by = c("Aberration" = "values")) %>%
      dplyr::mutate(SupportEv = paste0("TPM: ", tpm, ifelse(is.na(Pathway), "", paste0(", Pathway: ", Pathway)))) %>%
      filter(tpm < 10)  %>%
      dplyr::select(Aberration, Type,  Details, Variant_Properties, SupportEv, TargetValidation, Kids_First_Biospecimen_ID, sample_id, cohort, cohort_participant_id)
    
    # filter to amplifications that have TPM > 100
    myTableAmp <- myTable %>%
      filter(Type %in% c("Gain", "Amplification")) %>%
      left_join(rnaEvidence, by = c("Aberration" = "Gene")) %>%
      left_join(sigGeneSets, by = c("Aberration" = "values")) %>%
      dplyr::mutate(SupportEv = paste0("TPM: ", tpm, ifelse(is.na(Pathway), "", paste0(", Pathway: ", Pathway)))) %>%
      filter(tpm > 100)  %>%
      dplyr::select(Aberration, Type,  Details, Variant_Properties, SupportEv, TargetValidation, Kids_First_Biospecimen_ID, sample_id, cohort, cohort_participant_id)
    
    # mutation filters
    myTableMut <- myTable %>%
      filter(Type == "Mutation")  %>%
      left_join(rnaEvidence, by = c("Aberration" = "Gene")) %>%
      left_join(sigGeneSets, by = c("Aberration" = "values")) %>%
      dplyr::mutate(SupportEv = paste0("TPM: ", tpm, ifelse(is.na(Pathway), "", paste0(", Pathway: ", Pathway)))) %>%
      dplyr::select(Aberration, Type,  Details, Variant_Properties, SupportEv, TargetValidation, Kids_First_Biospecimen_ID, sample_id, cohort, cohort_participant_id)
    
    # fusion filters
    myTableFus <- myTable %>%
      filter(Type == "Fusion")
    if(nrow(myTableFus) > 0){
      myTableFus <- myTableFus %>%
        dplyr::mutate(Gene1 = sapply(Aberration, FUN = function(x) strsplit(x, "_")[[1]][[1]]),
                      Gene2 = sapply(Aberration, FUN = function(x) strsplit(x, "_")[[1]][[2]])) %>%
        left_join(rnaEvidence, by = c("Gene1" = "Gene"))  %>%
        dplyr::rename("Gene1_TPM" = tpm) %>%
        left_join(rnaEvidence, by = c("Gene2" = "Gene"))  %>%
        dplyr::rename("Gene2_TPM" = tpm)  %>%
        left_join(sigGeneSets, by = c("Gene1" = "values")) %>%
        left_join(sigGeneSets, by = c("Gene2" = "values")) %>%
        dplyr::mutate(SupportEv = paste0("Gene1_TPM: ", Gene1_TPM, ", ",
                                         "Gene2_TPM: ", Gene2_TPM,
                                         ifelse(is.na(Pathway.x), "", paste0(", Gene1_Pathway: ", Pathway.x)),
                                         ifelse(is.na(Pathway.y), "", paste0(", Gene2_Pathway: ", Pathway.y)))) %>%
        dplyr::select(Aberration, Type,  Details, Variant_Properties, SupportEv, TargetValidation, Kids_First_Biospecimen_ID, sample_id, cohort, cohort_participant_id)
    } else {
      myTableFus <- data.frame()
    }
    
    # merge all data
    myTable <- rbind(myTableAmp, myTableDel, myTableMut, myTableFus)
    myTable <- myTable %>%
      dplyr::rename("Supporting_Evidence" = "SupportEv") %>%
      unique()
  }
  
  return(myTable)
}
