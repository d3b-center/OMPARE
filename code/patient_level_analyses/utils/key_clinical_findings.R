# high confidence findings

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

key_clinical_findings <- function(snv_pattern, all_findings_output) {

  # use all findings as base
  myTable <- all_findings_output

  # expression is critical
  if(exists('rnaseq_analysis_output')){
    rnaEvidence <- rnaseq_analysis_output$TPM %>%
      rownames_to_column("Gene")
    
    # Get only significant sets with adj. pvalue < 0.05
    sigGeneSets <- rnaseq_analysis_output$pathways %>%
      filter(ADJ_P_VAL < 0.05)  %>%
      dplyr::select(Pathway, Direction) %>%
      inner_join(geneSetTS, by = c("Pathway" = "ind"))  %>%
      mutate(Pathway = paste0(Pathway, "(", Direction, ")")) %>%
      dplyr::select(Pathway, values) %>%
      group_by(values) %>%
      summarise(Pathway = toString(Pathway))
    
    # Supporting Evidence for Deletions
    myTableDel <- myTable %>%
      filter(Type == "Deletion") %>%
      left_join(rnaEvidence, by = c("Aberration" = "Gene")) %>%
      left_join(sigGeneSets, by = c("Aberration" = "values")) %>%
      mutate(SupportEv = paste0("TPM: ", !!as.name(sampleInfo$subjectID), ifelse(is.na(Pathway), "", paste0(", Pathway: ", Pathway)))) %>%
      filter(!!as.name(sampleInfo$subjectID) < 10)  %>%
      dplyr::select(Aberration, Type,  Details, Variant_Properties, SupportEv)
    
    # Supporting Evidence for Amplifications
    myTableAmp <- myTable %>%
      filter(Type == "Amplification") %>%
      left_join(rnaEvidence, by = c("Aberration" = "Gene")) %>%
      left_join(sigGeneSets, by = c("Aberration" = "values")) %>%
      mutate(SupportEv = paste0("TPM: ", !!as.name(sampleInfo$subjectID), ifelse(is.na(Pathway), "", paste0(", Pathway: ", Pathway)))) %>%
      filter(!!as.name(sampleInfo$subjectID) > 100)  %>%
      dplyr::select(Aberration, Type,  Details, Variant_Properties, SupportEv)
    
    # Supporting Evidence for Mutations Oncogene - Expression is listed, Pathway is significant
    myTableMut <- myTable %>%
      filter(Type == "Mutation")  %>%
      left_join(rnaEvidence, by = c("Aberration" = "Gene")) %>%
      left_join(sigGeneSets, by = c("Aberration" = "values")) %>%
      mutate(SupportEv = paste0("TPM: ", !!as.name(sampleInfo$subjectID), ifelse(is.na(Pathway), "", paste0(", Pathway: ", Pathway)))) %>%
      dplyr::select(Aberration, Type,  Details, Variant_Properties, SupportEv)
    
    # Supporting Evidence for Mutations Oncogene - Expression is listed, Pathway is significant
    myTableFus <- myTable %>%
      filter(Type == "Fusion") %>%
      mutate(Gene1 = sapply(Aberration, FUN = function(x) strsplit(x, "_")[[1]][[1]]),
             Gene2 = sapply(Aberration, FUN = function(x) strsplit(x, "_")[[1]][[2]])) %>%
      left_join(rnaEvidence, by = c("Gene1" = "Gene"))  %>%
      dplyr::rename("Gene1_TPM" = !!as.name(sampleInfo$subjectID)) %>%
      left_join(rnaEvidence, by = c("Gene2" = "Gene"))  %>%
      dplyr::rename("Gene2_TPM" = !!as.name(sampleInfo$subjectID))  %>%
      left_join(sigGeneSets, by = c("Gene1" = "values")) %>%
      left_join(sigGeneSets, by = c("Gene2" = "values")) %>%
      mutate(SupportEv = paste0("Gene1_TPM: ", Gene1_TPM, ", ",
                               "Gene2_TPM: ", Gene2_TPM,
                               ifelse(is.na(Pathway.x), "", paste0(", Gene1_Pathway: ", Pathway.x)),
                               ifelse(is.na(Pathway.y), "", paste0(", Gene2_Pathway: ", Pathway.y)))) %>%
      dplyr::select(Aberration, Type,  Details, Variant_Properties, SupportEv)
    
    # merge all data
    myTable <- rbind(myTableAmp, myTableDel, myTableMut, myTableFus)
    myTable <- myTable %>%
      dplyr::rename("Supporting Evidence" = "SupportEv") %>%
      unique()
  } else {
    myTable <- data.frame()
  }
  
  # add Ensembl ids and map to targetvalidation.org
  myTable <- myTable %>%
    left_join(expData %>% dplyr::select(gene_id, gene_symbol), by = c("Aberration" = "gene_symbol")) %>%
    mutate(TargetValidation = ifelse(is.na(gene_id), "", paste0('<a href = \"https://www.targetvalidation.org/target/',gene_id,'\">',gene_id,"</a>"))) %>%
    dplyr::select(-c(gene_id))
  
  return(myTable)
}
