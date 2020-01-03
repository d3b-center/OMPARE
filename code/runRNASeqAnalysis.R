##########################################
# Process RNA-Seq and Run Pathway Analysis
##########################################

# Function to return all results from RNA-Seq Analysis
suppressPackageStartupMessages(library(GSVA))
suppressPackageStartupMessages(library(GSEABase))

cancerGenes <- read.delim("data/Reference/CancerGeneList.tsv", stringsAsFactors = F)

runRNASeqAnalysis <- function(expData = NULL) {
  
  # Format gtex data ensembl ids
  gtexData$gene_id <- sapply(as.character(gtexData$gene_id), FUN=remDotStuff)
  rownames(gtexData) <- gtexData$gene_id
  gtexData <- gtexData[-1]

  # Merge GTEx and Patient data on common genes
  intGenesTmp <- intersect(rownames(gtexData), rownames(expData))
  mergeDF <- cbind(gtexData[intGenesTmp,], expData[intGenesTmp,"FPKM"])

  # Collapse to unique gene symbols
  # Matrix of GTEx and Patient FPKM data
  gtexGeneAnnot <- unique(gtexGeneAnnot[c("gene_id","gene_symbol")])
  gtexGeneAnnot$gene_id <- sapply(as.character(gtexGeneAnnot$gene_id), FUN=remDotStuff)
  rownames(gtexGeneAnnot) <- gtexGeneAnnot$gene_id
  mergeDF <- cbind(gtexGeneAnnot[rownames(mergeDF),2], mergeDF)
  mergeDF[,"meanExp"] <- rowMeans(mergeDF[-1]) 
  mergeDF <- mergeDF[order(-mergeDF[,"meanExp"]),]
  mergeDF <- mergeDF[!duplicated(mergeDF[,1]),]
  rownames(mergeDF) <- mergeDF[,1] 
  mergeDF <- mergeDF[-1]
  mergeDF <- mergeDF[-ncol(mergeDF)] 
  colnames(mergeDF)[ncol(mergeDF)] <- "SampleX"

  # Calculate Gene Outliers in Patient (top 20 Up and Down)
  getAllOutliers <- function(myMergeDF = mergeDF, getTop = 20, cancerGeneNames = cancerGenes$Gene) {
    
    # Filter in Patient: FPKM > 10 
    myMergeDF <- myMergeDF[myMergeDF$SampleX > 10,]
    
    # z-score and return only patient's value
    getZ <- function(x) {
      x <- log2(x+1)
      out <- (x-mean(x))/sd(x)
      return(out[length(out)])
    }
    output <- apply(myMergeDF, FUN=getZ, MARGIN=1)
    outputCanc <- output[intersect(names(output), cancerGeneNames)] # filter to cancer gene list
    outputDown <- sort(outputCanc)[1:getTop] # top 20 down
    outputUp <- sort(outputCanc, T)[1:getTop] # top 20 up
    
    outputUpDF <- data.frame(outputUp, myMergeDF[names(outputUp),"SampleX"])
    colnames(outputUpDF) <- c("Z_Score", "FPKM")
    outputDownDF <- data.frame(outputDown, myMergeDF[names(outputDown),"SampleX"])
    colnames(outputDownDF) <- c("Z_Score", "FPKM")
    
    return(list(output, rbind(outputUpDF, outputDownDF)))
  }
  geneAnalysisOut <- getAllOutliers()

  # Calculate Pathway Outliers-
  # Currently use Enrichment, but moving forward will use GSVA

  # Set Threshold
  thresh <- 1.5
  
  # Get Up and Down Genes
  tmpghj <- geneAnalysisOut[[1]]
  upGenes <- names(tmpghj)[which(geneAnalysisOut[[1]] > thresh)]
  downGenes <- names(tmpghj)[which(geneAnalysisOut[[1]] < (-1*thresh))]
  
  # If not enough genes take top 1000 
  if(length(upGenes) < 1000) {
    upGenes <- names(sort(geneAnalysisOut[[1]], decreasing = TRUE))[1:1000]
  }
  if(length(downGenes) < 1000) {
    downGenes <- names(sort(geneAnalysisOut[[1]], decreasing = FALSE))[1:1000]
  }
  
  # Code to run pathway analysis
  runHypGeom <- function(set, genes, n = 20000, universe = NULL) {
    
    if(!is.null(universe)) {
      set <- intersect(set, universe)
    }
    # number of white balls
    x <- length(intersect(genes, set))
    
    # white balls
    m <- length(genes)
    
    # black balls
    n2 <- n-m 
    
    # balls drawn from the urn 
    k <- length(set)
    
    out <- phyper(x-1, m, n2, k, lower.tail=F)
    setSize <- k
    overLap <- x
    numGenes <- m
    
    myRet <- c(setSize, numGenes, overLap, out) 
    return(myRet)
  }
  
  # Accessory for functional enrichment
  funcEnrichment <- function(genes, sets, qval=.25, numRet=5, myN=20000, myUniverse=NULL) {
    
    out <- lapply(sets, FUN = runHypGeom, genes = genes, n=myN, universe=myUniverse)
    out <- data.frame(out)
    out <- data.frame(t(out))
    out$ADJ_P_VAL <- p.adjust(out[,4], method="BH")
    colnames(out)[1:5] <- c("SET_SIZE", "NUM_GENES_INPUT", "OVERLAP", "P_VAL", "ADJ_P_VALUE")
    return(out)
    
  }
  
  upPathways <- funcEnrichment(upGenes, hallMarkSets, qval=1, myN=25000, myUniverse=rownames(mergeDF))
  upPathways <- upPathways[order(upPathways$P_VAL),]
  upPathways[,"Direction"] <- "Up"
  upPathways[,"Pathway"] <- rownames(upPathways)
  
  downPathways <- funcEnrichment(downGenes, hallMarkSets, qval=1, myN=25000, myUniverse=rownames(mergeDF))
  downPathways <- downPathways[order(downPathways$P_VAL),]
  downPathways[,"Direction"] <- "Down"
  downPathways[,"Pathway"] <- rownames(downPathways)
  
  pathwayAnalysisOut <- list(list("UpGenes" = upGenes, "DownGenes" = downGenes), rbind(upPathways, downPathways))

  # Final Output
  finalOut <- list()
  finalOut$geneAnalysis <- geneAnalysisOut
  finalOut$pathwayAnalysis <- pathwayAnalysisOut
  finalOut$FPKM <- mergeDF["SampleX"]
  return(finalOut)
}
