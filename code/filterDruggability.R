#####################################
# Purpose: Pull Annotation from dgiDB
#####################################

filterDruggability <- function(data = dgidb){
  
  # Filter to known interaction types
  dataOut <- data %>%
    filter(interaction_types != "" & drug_name != "" & gene_name != "") %>%
    dplyr::select(gene_name, drug_name) %>%
    unique() %>%
    group_by(gene_name) %>% 
    dplyr::summarize(Drugs=paste(drug_name, collapse=", ")) %>%
    as.data.frame()
    
  return(dataOut)
}
