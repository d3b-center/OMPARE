cnv_plot <- function(myCnvData = cnvRatioData, fname) {
  myCnvData <- myCnvData[,c('Chromosome', 'Start', 'CopyNumber')]
  myCnvData.seg <- pcf(data = myCnvData, verbose = FALSE)
  png(filename = fname, height = 4, width = 10, units = "in", res = 300)
  plotGenome(data = myCnvData, segments = myCnvData.seg)
  dev.off()
}