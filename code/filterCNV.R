########################
# Function to filter CNV
########################

cancerGenes <- read.delim("data/Reference/CancerGeneList.tsv", stringsAsFactors = F)

filterCNV <- function(myCNVData = cnvGenes, myCancerGenes = cancerGenes, myTSGenes = tsgGenes, cutoffHigh = 4, cutoffLow = 1) {
  cnvDataFilt <- myCNVData
  
  # Filter by Oncogenes/TSGs as well as copy number cutoff
  myTSGenes <- as.character(myTSGenes$GeneSymbol)
  myOncogenes <- setdiff(as.character(myCancerGenes$Gene), myTSGenes)
  cnvDataFiltUp <- cnvDataFilt %>%
    filter(CNA > cutoffHigh & Gene %in% myOncogenes)
  cnvDataFiltDown <- cnvDataFilt %>%
    filter(CNA < cutoffLow & Gene %in% myTSGenes)
  cnvDataFilt <- rbind(cnvDataFiltUp, cnvDataFiltDown)
  
  return(cnvDataFilt)
}
