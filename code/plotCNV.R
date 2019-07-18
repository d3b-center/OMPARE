#######################
# Function for CNV View
#######################

plotCNV <- function(myCnvData = cnvData) {
  colnames(myCnvData) <- c("Chr", "StartBP", "EndBP", "AbsCNV", "Type")
  myCnvData[,"log2CNV"] <- log2((myCnvData[,"AbsCNV"]+1))
  myCnvData[,"MedianBP"] <- round((myCnvData[,"StartBP"]+myCnvData[,"EndBP"])/2)
  myCnvData <- myCnvData[,c("Chr", "MedianBP", "log2CNV")]
  single.seg <- pcf(data=myCnvData,gamma=12,verbose=FALSE)
  plotGenome(myCnvData, single.seg)
}
