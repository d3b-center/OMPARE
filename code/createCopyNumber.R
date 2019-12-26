#######################################
# Copy number by gene instead of region
#######################################

createCopyNumber <- function(){
  
  # Chromosome Map
  chrMap <- chrMap[!chrMap[,1]=="",]
  getCNV <- function(x) {
    tmpChr <- x[[1]]
    tmpStart <- as.numeric(trimws(x[[2]]));
    tmpEnd <- as.numeric(trimws(x[[3]]));
    tmpValue <- as.numeric(trimws(x[[4]]));
    tmpCNA <- chrMap[chrMap[,"Chromosome.scaffold.name"] == tmpChr,];
    tmpCNA <- tmpCNA[((tmpCNA[,2] > tmpStart) & (tmpCNA[,2] < tmpEnd)),1]
    if(length(tmpCNA)>0) {
      tmpCNA <- data.frame(tmpCNA, tmpValue);
    }
    return(tmpCNA);
  }
  
  outList <- apply(cnvData, FUN=getCNV, MARGIN=1)
  output <- do.call("rbind", outList);
  missingGenes <- setdiff(chrMap[,1], output[,1]);
  missingGenes <- data.frame(missingGenes, 2)
  colnames(missingGenes) <- colnames(output)
  output <- rbind(output, missingGenes);
  colnames(output) <- c("Gene", "CNA")
  return(output);
}










