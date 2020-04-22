# Disease Specific Information
# Gets information from allFindingsTable and removes Pathway information
diseaseSpecificInformation <- function() {
  tmpGeneFindings <- allFindingsTable()
  if(nrow(tmpGeneFindings) > 0){
    tmpGeneFindings <- tmpGeneFindings[!grepl("Pathway", tmpGeneFindings[,"Type"]),]
    diseaseSpecificFields
    
    # Check everything
    getStatus <- function(x) {
      tmpGenes <- as.character(x[[3]]) # Gene
      tmpGenes <- trimws(strsplit(tmpGenes, ",")[[1]]) # For multiple genes, split into vector
      tmpOut <- sapply(paste0('^', tmpGenes), FUN = grepl, x = tmpGeneFindings[,1]) # find genes in all findings table (without pathway info)
      paste(paste(tmpGeneFindings[as.logical(rowSums(tmpOut)),1], ":", tmpGeneFindings[as.logical(rowSums(tmpOut)),2], sep=""), collapse=", ")  # get Aberration and Type info 
    }
    
    diseaseSpecificFields[,"Value"] <- apply(diseaseSpecificFields, FUN=getStatus, MARGIN=1)
    diseaseSpecificFields[diseaseSpecificFields[,"Value"]==":","Value"] <- "Normal" # if no information is found from all findings, just assign Normal
    diseaseSpecificFields <- diseaseSpecificFields[,c("Field_name", "Value")]
  } else {
    diseaseSpecificFields <- data.frame()
  }
  return(diseaseSpecificFields)
}
