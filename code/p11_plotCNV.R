#######################
# Function for CNV View
#######################

plotCNV <- function(myCnvData = cnvRatioData) {
  myCnvData <- myCnvData[,c('Chromosome','Start','CopyNumber')]
  myCnvData.seg <- pcf(data=myCnvData, verbose=FALSE)
  plotGenome(myCnvData, myCnvData.seg)
}