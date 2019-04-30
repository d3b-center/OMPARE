###############################
#Purpose: Copy number by gene instead of region
#Date: 2/11/2019
#Author: Pichai Raman
###############################



createCopyNumber <- function()
{


#Read copy number data
cnvData <- read.delim("../data/CBTTC-HGG/CNV/5c3eace5-950a-4a05-81ed-5c04b4a0a367.CNVs", header=F);

#Chromosome Map
chrMap <- read.delim("../data/Reference/mart_export_genechr_mapping.txt", stringsAsFactors =F)
chrMap <- chrMap[!chrMap[,1]=="",]
getCNV <- function(x)
{
  tmpChr <- x[[1]]
  tmpStart <- as.numeric(trimws(x[[2]]));
  tmpEnd <- as.numeric(trimws(x[[3]]));
  tmpValue <- as.numeric(trimws(x[[4]]));
  tmpCNA <- chrMap[chrMap[,"Chromosome.scaffold.name"]==tmpChr,];
  tmpCNA <- tmpCNA[((tmpCNA[,2]>tmpStart)&(tmpCNA[,2]<tmpEnd)),1]
  if(length(tmpCNA)>0)
  {
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










