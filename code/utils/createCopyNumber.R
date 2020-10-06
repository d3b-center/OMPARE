#######################################
# Copy number by gene instead of region
#######################################

createCopyNumber <- function(cnvData, ploidy = NULL){
  
  # Chromosome Map
  chrMap <- chrMap %>%
    filter(hgnc_symbol != "")
  
  getCNV <- function(x) {
    tmpChr <- x['chr']
    tmpStart <- as.numeric(x['start'])
    tmpEnd <- as.numeric(x['end'])
    tmpValue <- as.numeric(x['copy.number'])
    tmpStatus <- as.character(x['status'])
    tmpPval <- scientific(as.numeric(x['WilcoxonRankSumTestPvalue']), digits = 2) 
    tmpCNA <- chrMap[chrMap$chromosome == tmpChr,]
    tmpCNA <- tmpCNA[((tmpCNA$gene_start > tmpStart) & (tmpCNA$gene_end < tmpEnd)), 'hgnc_symbol']
    if(length(tmpCNA) > 0) {
      tmpCNA <- data.frame(tmpCNA, tmpValue, tmpStatus, tmpPval)
    }
    return(tmpCNA)
  }
  
  outList <- apply(cnvData, FUN = getCNV, MARGIN=1)
  output <- do.call("rbind", outList)
  # if ploidy is not supplied explicitly, use value corresponding to neutral status
  if(is.null(ploidy)){
    ploidy <- unique(output[which(output$tmpStatus == "neutral"),'tmpValue'])
  }
  if(length(ploidy) > 1){
    ploidy <- min(ploidy)
  }
  # if there is no CN entry for neutral status, then do this
  if(length(ploidy) == 0){
    ploidy <- min(output$tmpValue[output$tmpStatus == "gain"])-1
  }
  missingGenes <- setdiff(chrMap$hgnc_symbol, output$tmpCNA)
  missingGenes <- data.frame(missingGenes, ploidy, 'neutral', 1)
  colnames(missingGenes) <- colnames(output)
  output <- rbind(output, missingGenes)
  colnames(output) <- c("Gene", "CNA", "Status", "Pvalue")
  output <- unique(output)
  return(output)
}










