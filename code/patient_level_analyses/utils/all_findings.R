# all findings table 

all_findings <- function(snv_caller) {
  # Add Mutations
  if(exists('mutDataAnnot')){
    # annotate mutations, add Aberration and Details
    tmpMut <- mutDataAnnot %>%
      mutate(Aberration = Hugo_Symbol) %>%
      filter(HGVSp_Short != "") %>%
      mutate(Details = paste0('Variant: ', Variant_Classification, " | HGVSp: ", HGVSp_Short)) 
    
    # only keep Pfam and SMART domains, collapse DOMAINS to string
    tmpMut <- tmpMut %>%
      mutate(DOMAINS = ifelse(DOMAINS == "", NA, strsplit(as.character(DOMAINS), ","))) %>% 
      unnest(DOMAINS) %>%
      mutate(match = str_detect(DOMAINS, 'Pfam_domain|SMART_domains')) %>%
      mutate(DOMAINS = ifelse(match, DOMAINS, NA)) %>%
      group_by(Aberration, Type, Details) %>%
      mutate(DOMAINS = toString(na.omit(DOMAINS))) %>%
      unique()
    
    # add COSMIC
    tmpMut <- tmpMut %>%
      mutate(Existing_variation = str_detect(Existing_variation, 'COSM')) %>%
      mutate(Existing_variation = ifelse(Existing_variation, "Cosmic_Variant", ""))
    
    # now add Variant Properties depending on snv_caller
    if(snv_caller == "all"){
      tmpMut <- tmpMut %>%
        dplyr::select(Aberration, Type, Details, SIFT, DOMAINS, Existing_variation) %>%
        unique()
    } else {
      tmpMut <- tmpMut %>%
        dplyr::select(Aberration, Type, Details, t_alt_count, t_depth, SIFT, DOMAINS, Existing_variation) %>%
        unique()
    }
    
    # collapse into Variant_Properties
    tmpMut <- reshape2::melt(as.data.frame(tmpMut), id.vars = c('Aberration', 'Type', 'Details'))
    tmpMut <- tmpMut %>%
      group_by(Aberration, Type, Details, variable) %>%
      mutate(Variant_Properties = paste0(variable,":", value))  %>%
      group_by(Aberration, Type, Details) %>%
      summarise(Variant_Properties = paste(Variant_Properties, collapse = '; ')) %>%
      unique()
  } else {
    tmpMut <- data.frame()
  }
  
  # Add Fusions
  if(exists('fusData')){
    tmpFus <- fusData %>%
      mutate(Aberration = X.fusion_name,
             Type = "Fusion",
             Details = paste0("Fusion Type: ", Splice_type),
             Variant_Properties = "") %>%
      dplyr::select(Aberration, Type, Details, Variant_Properties)
  } else {
    tmpFus <- data.frame()
  }
  
  # Add Copy Number
  if(exists('cnvDataFilt')){
    tmpCnv <- cnvDataFilt %>%
      rowwise() %>%
      mutate(Aberration = hgnc_symbol,
             Type = status,
             Details = paste0("Copy Number: ", copy.number, 
                              " | Pos: ", chr, ":", start, "-", end),
             Variant_Properties = NA) %>%
      dplyr::select(Aberration, Type, Details, Variant_Properties)
  } else {
    tmpCnv <- data.frame()
  }
  
  # Add Expression
  if(exists('expData')){
    tmpExp <- rnaseq_analysis_output$diffexpr.top20 %>%
      rownames_to_column("Aberration") %>%
      mutate(Type = c(rep("Outlier-High (mRNA)", 20), rep("Outlier-Low (mRNA)", 20)),
             Details = paste0("logFC: ", round(logfc, 2)," | TPM: ", tpm),
             Variant_Properties = "")  %>%
      dplyr::select(Aberration, Type, Details, Variant_Properties)
  } else {
    tmpExp <- data.frame()
  }
  
  # Add Pathway
  if(exists('expData')){
    tmpPath <- rnaseq_analysis_output$pathways
    tmpPathUp <- tmpPath %>%
      filter(direction == "up") %>%
      mutate(Type = "Pathway Up") %>%
      arrange(padj) %>%
      slice_head(n = 20)
    tmpPathDown <- tmpPath %>%
      filter(direction == "down") %>%
      mutate(Type = "Pathway Down") %>%
      arrange(padj) %>%
      slice_head(n = 20)
    tmpPath <- rbind(tmpPathUp, tmpPathDown)
    tmpPath <- tmpPath %>%
      mutate(Aberration = pathway,
             Details = paste0("Adj. P-Value: ", formatC(padj, format = "e", digits = 2)),
             Variant_Properties = "") %>%
      dplyr::select(Aberration, Type, Details, Variant_Properties)
  } else {
    tmpPath <- data.frame()
  }
  
  # Merge Together
  allFindingsDF <- unique(rbind(tmpMut, tmpFus, tmpCnv, tmpExp, tmpPath))
  
  # add Ensembl ids and map to targetvalidation.org
  myTable <- allFindingsDF %>%
    left_join(expData %>% dplyr::select(gene_id, gene_symbol), by = c("Aberration" = "gene_symbol")) %>%
    mutate(TargetValidation = ifelse(is.na(gene_id), "", paste0('<a href = \"https://www.targetvalidation.org/target/',gene_id,'\">',gene_id,"</a>"))) %>%
    dplyr::select(-c(gene_id)) 
  
  return(myTable)
}
