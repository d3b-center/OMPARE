# filter fusions

filter_fusions <- function(fusion_data, myCancerGenes = cancer_genes, method, myJunctionReads = 2) {
  

  # format column names
  if(method == "star-fusion"){
    fusion_data <- fusion_data %>%
      mutate(gene1 = gsub("--.*", "", FusionName),
             gene2 = gsub(".*--", "", FusionName),
             fusion_name = gsub('--','_', FusionName),
             type = PROT_FUSION_TYPE)
  } else {
    fusion_data <- fusion_data %>%
      mutate(fusion_name = paste0(as.character(gene1), "_", as.character(gene2)))
  }
  
  # filter by number of reads (star) or confidence (arriba)
  if(method == "star-fusion"){
    fusion_data <- fusion_data %>% 
      filter(JunctionReadCount > myJunctionReads)
  } else {
    fusion_data <- fusion_data %>%
      filter(confidence != "low")
  }
  
  # filter by cancer genes (from annoFuse)
  myCancerGenes <- as.character(myCancerGenes$Gene_Symbol)
  fusion_data <- fusion_data %>%
    filter(gene1 %in% myCancerGenes | gene2 %in% myCancerGenes) %>%
    dplyr::select(fusion_name, gene1, gene2, annots, type, tumor_id) %>%
    unique()
  
  # fix missing values
  fusion_data <- fusion_data %>%
    mutate(annots = ifelse(annots == "[]", "", annots),
           type = ifelse(type == ".", "", type))
    
  return(fusion_data)
}
