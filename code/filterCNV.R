########################
# Function to filter CNV
########################

filterCNV <- function(myCNVData = cnvGenes, myCancerGenes = cancerGenes, myTSGenes = tsgGenes, cutoffHigh = 4, cutoffLow = 1) {
  cnvDataFilt <- myCNVData
  
  # Filter by Oncogenes/TSGs as well as copy number cutoff
  myTSGenes <- cancerGenes %>%
    filter(type == "TumorSuppressorGene") %>%
    .$Gene_Symbol %>%
    unique()
  myOncogenes <- cancerGenes %>%
    filter(type %in% c("Oncogene", "OncoKB")) %>%
    .$Gene_Symbol %>%
    unique()
  myOncogenes <- setdiff(myOncogenes, myTSGenes)
  cnvDataFiltUp <- cnvDataFilt %>%
    filter(CNA > cutoffHigh & Gene %in% myOncogenes)
  cnvDataFiltDown <- cnvDataFilt %>%
    filter(CNA < cutoffLow & Gene %in% myTSGenes)
  cnvDataFilt <- rbind(cnvDataFiltUp, cnvDataFiltDown)
  
  return(cnvDataFilt)
}
