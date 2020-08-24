# Disease Specific Information
# Gets information from allFindingsTable and removes Pathway information
diseaseSpecificInformation <- function() {
  
  tmpGeneFindings <- allFindingsTable() %>%
    filter(!Type %in% c("Pathway Up", "Pathway Down")) %>%
    mutate(Aberration = ifelse(Type == "Mutation", yes = paste0(Aberration, ": ", gsub('.*HGVSp: ','', Details)), no = Aberration))

  # Check everything
  getStatus <- function(x) {
    tmpGenes <- as.character(x[[3]]) # Gene
    tmpGenes <- trimws(strsplit(tmpGenes, ",")[[1]]) # For multiple genes, split into vector
    tmpOut <- sapply(paste0('^', tmpGenes), FUN = grepl, x = tmpGeneFindings[,1]) # find genes in all findings table (without pathway info)
    paste(paste(tmpGeneFindings[as.logical(rowSums(tmpOut)),1], " (", tmpGeneFindings[as.logical(rowSums(tmpOut)),2], ")", sep=""), collapse=", ")  # get Aberration and Type info 
  }
  
  diseaseSpecificFields <- diseaseSpecificFields %>% 
    mutate(Value = apply(diseaseSpecificFields, FUN=getStatus, MARGIN=1)) %>%
    mutate(Value = ifelse(Value == " ()", "Normal", Value)) %>%
    dplyr::select(Field_name, Value)
  
  return(diseaseSpecificFields)
}
