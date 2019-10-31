########################
# Function to filter CNV
########################

filterCNV <- function(myCNVData=cnvGenes, myCancerGenes=cancerGenes, myTSGenes=tsgGenes, cutoffHigh=4, cutoffLow=1) {
  cnvDataFilt <- myCNVData
  
  # Filter by Cancer Gene Census
  myTSGenes <- as.character(myTSGenes[,2])
  myOncogenes <- setdiff(as.character(myCancerGenes[,1]), myTSGenes)
  
  cnvDataFiltUp <- cnvDataFilt[cnvDataFilt[,2]>cutoffHigh,]
  cnvDataFiltUp <- cnvDataFiltUp[cnvDataFiltUp[,1]%in%myOncogenes,]
  
  cnvDataFiltDown <- cnvDataFilt[cnvDataFilt[,2]<cutoffLow,]
  cnvDataFiltDown <- cnvDataFiltDown[cnvDataFiltDown[,1]%in%myTSGenes,]
  cnvDataFilt <- rbind(cnvDataFiltUp, cnvDataFiltDown)
  
  return(cnvDataFilt)
  
}
