#######################
# Function for CNV View
#######################

# plotCNV <- function(myCnvData = cnvData) {
#   colnames(myCnvData) <- c("Chr", "StartBP", "EndBP", "AbsCNV", "Type")
#   # myCnvData[,"log2CNV"] <- log2((myCnvData[,"AbsCNV"]+1))
#   myCnvData[,"CNV"] <- myCnvData[,"AbsCNV"]
#   myCnvData[,"MedianBP"] <- round((myCnvData[,"StartBP"]+myCnvData[,"EndBP"])/2)
#   myCnvData <- myCnvData[,c("Chr", "MedianBP", "CNV")]
#   single.seg <- pcf(data=myCnvData, gamma=12, verbose=FALSE)
#   plotGenome(myCnvData, single.seg)
# }

plotCNV <- function(myCnvData = cnvRatioData) {
  myCnvData <- myCnvData[,c('Chromosome','Start','CopyNumber')]
  myCnvData.seg <- pcf(data=myCnvData, verbose=FALSE)
  plotGenome(myCnvData, myCnvData.seg)
}