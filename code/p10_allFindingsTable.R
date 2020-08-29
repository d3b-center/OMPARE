#####################
# All Findings Table
#####################

source('code/filterCNV.R')
source('code/annotateMutations.R')

allFindingsTable <- function(snv_pattern) {
  # Add Mutations
  if(exists('mutData')){
    # annotate mutations, add Aberration and Details
    tmpMut <- annotateMutations() %>%
      mutate(Aberration = Hugo_Symbol) %>%
      filter(HGVSp_Short != "") %>%
      mutate(Details = paste0('Variant: ', Variant_Classification, " | HGVSp: ", HGVSp_Short)) 
    
    # remove all PDB-ENSP from DOMAINS, collapse DOMAINS to string
    tmpMut <- tmpMut %>%
      mutate(DOMAINS = ifelse(DOMAINS == "", NA, strsplit(as.character(DOMAINS), ","))) %>% 
      unnest(DOMAINS) %>%
      mutate(match = str_detect(DOMAINS, 'PDB-ENSP')) %>%
      mutate(DOMAINS = ifelse(match, NA, DOMAINS)) %>%
      group_by(Aberration, Type, Details) %>%
      mutate(DOMAINS = toString(na.omit(DOMAINS))) %>%
      unique()
    
    # now add Variant Properties depending on snv_pattern
    if(snv_pattern == "all"){
      tmpMut <- tmpMut %>%
        dplyr::select(Aberration, Type, Details, SIFT, DOMAINS) %>%
        unique()
    } else {
      tmpMut <- tmpMut %>%
        dplyr::select(Aberration, Type, Details, t_alt_count, t_depth, SIFT, DOMAINS) %>%
        unique()
    }
    
    # collapse into Variant_Properties
    tmpMut <- melt(tmpMut, id.vars = c('Aberration', 'Type', 'Details'))
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
  if(exists('cnvGenes')){
    tmpCnv <- filterCNV() %>%
      mutate(Aberration = Gene,
             Type = ifelse(Status == "gain", "Amplification", "Deletion"),
             Details = paste0("Copy Number Value: ", CNA),
             Variant_Properties = "") %>%
      dplyr::select(Aberration, Type, Details, Variant_Properties)
  } else {
    tmpCnv <- data.frame()
  }
  
  # Add Expression
  if(exists('expData')){
    tmpExp <- RNASeqAnalysisOut[[1]][[2]] %>%
      rownames_to_column("Aberration") %>%
      mutate(Type = c(rep("Outlier-High (mRNA)", 20), rep("Outlier-Low (mRNA)", 20)),
             Details = paste0("Z-score: ", round(Z_Score, 2)," | TPM: ",TPM),
             Variant_Properties = "")  %>%
      dplyr::select(Aberration, Type, Details, Variant_Properties)
  } else {
    tmpExp <- data.frame()
  }
  
  # Add Pathway
  if(exists('expData')){
    tmpPath <- RNASeqAnalysisOut[[2]][[2]]
    tmpPathUp <- tmpPath %>%
      filter(Direction == "Up") %>%
      arrange(ADJ_P_VALUE) %>%
      top_n(20)
    tmpPathDown <- tmpPath %>%
      filter(Direction == "Down") %>%
      arrange(ADJ_P_VALUE) %>%
      top_n(20)
    tmpPath <- rbind(tmpPathUp, tmpPathDown)
    tmpPath <- tmpPath %>%
      mutate(Aberration = Pathway,
             Type = c(rep("Pathway Up", 20), rep("Pathway Down", 20)),
             Details = paste0("P-Value: ", formatC(P_VAL, format = "e", digits = 2)),
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
  
  return(allFindingsDF)
}
