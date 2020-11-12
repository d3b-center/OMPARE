cnv_plot <- function(myCnvData = cnvGenes, fname) {
  myCnvData <- myCnvData %>%
    dplyr::select(chr, start, copy.number) %>%
    as.data.frame()
  myCnvData.seg <- copynumber::pcf(data = myCnvData, verbose = FALSE)
  png(filename = fname, height = 4, width = 10, units = "in", res = 300)
  plotGenome(data = myCnvData, segments = myCnvData.seg)
  dev.off()
}
