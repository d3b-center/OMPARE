# gets information from all_findings_table and removes pathway information

suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
})

disease_specific_information <- function(disease_specific_fields, all_findings_output) {
  
  # only up/down genes, fusions, mutations, vus, cnv
  tmpGeneFindings <- all_findings_output %>%
    filter(!Type %in% c("Pathway Up", "Pathway Down")) %>%
    mutate(Aberration = ifelse(Type == "Mutation", yes = paste0(Aberration, ": ", gsub('.*HGVSp: ','', Details)), no = Aberration)) %>%
    as.data.frame()
  
  # Check everything
  getStatus <- function(x, tmpGeneFindings) {
    tmpGenes <- as.character(x[['Gene']]) # Gene
    tmpGenes <- trimws(strsplit(tmpGenes, ",")[[1]]) # For multiple genes, split into vector
    tmpOut <- sapply(paste0('^', tmpGenes), FUN = grepl, x = tmpGeneFindings$Aberration) # find genes in all findings table (without pathway info)
    paste(paste(tmpGeneFindings[as.logical(rowSums(tmpOut)),'Aberration'], " (", 
                tmpGeneFindings[as.logical(rowSums(tmpOut)), 'Type'], ")", sep=""), 
          collapse=", ")  # get Aberration and Type info 
  }
  
  disease_specific_fields <- disease_specific_fields %>% 
    mutate(Value = apply(disease_specific_fields, FUN = function(x) getStatus(x = x, tmpGeneFindings = tmpGeneFindings), MARGIN=1)) %>%
    mutate(Value = ifelse(Value == " ()", "Normal", Value)) %>%
    dplyr::select(Field_name, Value)
  
  return(disease_specific_fields)
}
