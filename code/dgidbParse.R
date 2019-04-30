###########################################
#Purpose: Pull Annotation from dgiDB
#Author: Pichai Raman
#Date: 3/21/2019
###########################################


filterDruggability <- function()
{
  
  library("tidyverse")
  data <- read.delim("../data/Reference/DGIdb.txt", stringsAsFactors = F)
  
  #Filter to known interaction types
  data <- data[data[,"interaction_types"]!="",]
  data <- data[data[,"drug_name"]!="",]
  data <- unique(data[,c("gene_name", "drug_name")]);
  dataOut <- data %>% group_by(gene_name) %>% summarize(Drugs=paste(drug_name, collapse=", "))
  dataOut <- data.frame(dataOut);
  return(dataOut);
  
}