# Author: Komal S. Rathi
# Accessory functions for pathway enrichment

# z-score and return only patient's value (i.e. last column)
get_zscore <- function(x) {
  x <- log2(x+1)
  out <- (x-mean(x))/sd(x)
  return(out[length(out)])
}

# get up/down genes using z-score
get_all_outliers <- function(myMergeDF, getTop, thresh, comparison, cancerGeneNames) {
  
  # Filter in Patient: TPM > 10
  myMergeDF <- myMergeDF[myMergeDF$SampleX > 10,] 
  
  # z-score
  output <- apply(myMergeDF, FUN = get_zscore, MARGIN = 1) 
  
  # full data
  # combine TPM and z-score for sample of interest
  genes.df <- data.frame(Z_Score = output, TPM = myMergeDF[names(output),"SampleX"])
  genes.df$DE <- ""
  genes.df$DE[which(genes.df$Z_Score < (-1*thresh))] <- "Down"
  genes.df$DE[which(genes.df$Z_Score > thresh)] <- "Up"
  genes.df$Comparison <- comparison
  genes.df$CancerGene <- ifelse(rownames(genes.df) %in% cancerGeneNames, TRUE, FALSE)
  
  # top 20 Up/Down genes 
  outputCanc <- output[intersect(names(output), cancerGeneNames)] # filter to cancer gene list
  outputDown <- sort(outputCanc)[1:getTop] # top 20 down
  outputUp <- sort(outputCanc, T)[1:getTop] # top 20 up
  outputUpDF <- data.frame(Z_Score = outputUp, TPM = myMergeDF[names(outputUp),"SampleX"])
  outputDownDF <- data.frame(Z_Score = outputDown, TPM = myMergeDF[names(outputDown),"SampleX"])
  diffexpr.top20 = rbind(outputUpDF, outputDownDF)
  
  return(list(expr.genes.z.score = output,
              diffexpr.top20 = diffexpr.top20,
              genes.df = genes.df))
}

# get up/down pathways using hypergeometric test
run_hypgeom_test <- function(set, genes, n = 20000, universe = NULL) {
  
  if(!is.null(universe)) {
    set <- intersect(set, universe)
  }
  
  # number of white balls
  x <- length(intersect(genes, set))
  ngenes <- as.character(toString(intersect(genes, set)))
  m <- length(genes)
  
  # black balls
  n2 <- n-m
  
  # balls drawn from the urn
  k <- length(set)
  
  out <- phyper(x-1, m, n2, k, lower.tail=F)
  setSize <- k
  overLap <- x
  numGenes <- m
  
  myRet <- c(setSize, numGenes, overLap, out, ngenes)
  return(myRet)
}

# functional enrichment
func_enrichment <- function(genes, sets, qval = 0.25, numRet = 5, myN = 20000, myUniverse = NULL) {
  
  out <- lapply(sets, FUN = run_hypgeom_test, genes = genes, n = myN, universe = myUniverse)
  out <- data.frame(out)
  out <- data.frame(t(out))
  out[,4] <- as.numeric(as.character(out[,4]))
  out$ADJ_P_VAL <- p.adjust(out[,4], method="BH")
  colnames(out)[1:5] <- c("SET_SIZE", "NUM_GENES_INPUT", "OVERLAP", "P_VAL", "Gene_name")
  out$Pathway <- rownames(out)
  return(out)
  
}

# function to tabulate DE and Pathway results
run_rnaseq_analysis <- function(exp.data, refData, thresh = 2, comparison, single_sample = FALSE, sample_name = "SampleX") {
  
  # current sample of interest
  smp <- unique(as.character(exp.data$Sample))
  
  # add genes as rownames and remove sample name
  exp.data <- exp.data %>%
    remove_rownames() %>%
    column_to_rownames('Gene') %>%
    dplyr::select(-c(Sample))
  
  # now remove the current sample from refData (required when we are doing sample from the same cohort)
  refData <- refData[,grep(smp, colnames(refData), invert = T, value = T)]
  
  # Merge refData and sample of interest data on common genes
  mergeDF <- refData %>% rownames_to_column('gene') %>% 
    inner_join(exp.data %>% 
                 rownames_to_column('gene'), by = 'gene') %>% 
    column_to_rownames('gene')
  colnames(mergeDF)[ncol(mergeDF)] <- "SampleX" # sample of interest
  
  # get gene outliers (top 20 Up and Down)
  geneAnalysisOut <- get_all_outliers(myMergeDF = mergeDF, getTop = 20, thresh = thresh, comparison = comparison, cancerGeneNames = cancerGenes$Gene_Symbol)
  
  # get up/down genes (using thresh as z-score cutoff)
  expr.genes.z.score <- geneAnalysisOut$expr.genes.z.score
  upGenes <- names(expr.genes.z.score)[which(expr.genes.z.score > thresh)]
  downGenes <- names(expr.genes.z.score)[which(expr.genes.z.score < (-1*thresh))]
  diff.genes = list("UpGenes" = upGenes, "DownGenes" = downGenes)
  
  # up pathways (using adj. pvalue < 0.05 cutoff)
  upPathways <- func_enrichment(upGenes, geneSet, qval = 1, myN = 25000, myUniverse = rownames(mergeDF))
  upPathways <- upPathways %>%
    mutate(Direction = "Up") %>%
    filter(ADJ_P_VAL < 0.05) %>%
    arrange(ADJ_P_VAL)
  
  # down pathways (using adj. pvalue < 0.05 cutoff)
  downPathways <- func_enrichment(downGenes, geneSet, qval = 1, myN = 25000, myUniverse = rownames(mergeDF))
  downPathways <- downPathways %>%
    mutate(Direction = "Down") %>%
    filter(ADJ_P_VAL < 0.05) %>%
    arrange(ADJ_P_VAL)
  
  # full pathway dataframe
  pathway.df <- upPathways %>% 
    rbind(downPathways) %>% 
    dplyr::select(Pathway, Gene_name, SET_SIZE, NUM_GENES_INPUT, OVERLAP, P_VAL, ADJ_P_VAL, Direction) %>% 
    mutate(Comparison = comparison)
  
  # expand gene names and get a map of genes + pathways
  pathway.df.exp <- pathway.df %>%
    separate_rows(Gene_name, convert = TRUE) %>%
    dplyr::select(Pathway, Gene_name)
  
  # create empty df just so you can map to gene.df below
  if(nrow(pathway.df.exp) == 0){
    pathway.df.exp <- data.frame(Pathway = "", Gene_name = "")
  }
  
  # map pathways to individual genes
  genes.df <- geneAnalysisOut$genes.df  %>%
    rownames_to_column("Gene_name") %>%
    left_join(pathway.df.exp, by = "Gene_name") %>%
    filter(DE != "") %>%
    group_by(Gene_name) %>%
    mutate(Pathway = toString(Pathway)) %>%
    unique() 

  if(single_sample == TRUE){
    colnames(mergeDF)[ncol(mergeDF)] <- sample_name
    finalOut <- list(pathways = pathway.df, 
                     genes = genes.df,
                     diff.genes = diff.genes,
                     TPM = mergeDF[ncol(mergeDF)],
                     expr.genes.z.score = geneAnalysisOut$expr.genes.z.score,
                     diffexpr.top20 = geneAnalysisOut$diffexpr.top20)
  } else {
    finalOut <- list(pathways = pathway.df, 
                     genes = genes.df)
  }
  
  return(finalOut)
}
