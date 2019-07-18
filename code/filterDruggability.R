#####################################
# Purpose: Pull Annotation from dgiDB
#####################################

library("tidyverse")

filterDruggability <- function(data = dgidb) {
  
  # Filter to known interaction types
  data <- data[data[,"interaction_types"]!="",]
  data <- data[data[,"drug_name"]!="",]
  data <- unique(data[,c("gene_name", "drug_name")]);
  dataOut <- data %>% group_by(gene_name) %>% summarize(Drugs=paste(drug_name, collapse=", "))
  dataOut <- data.frame(dataOut);
  return(dataOut);
  
}