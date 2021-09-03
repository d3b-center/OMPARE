# Author: Komal S. Rathi
# fgsea based pathway enrichment

suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(fgsea)
})

# z-score and return only patient's value (i.e. last column)
get_zscore <- function(x) {
  x <- log2(x+1)
  out <- (x-mean(x))/sd(x)
  return(out[length(out)])
}

# get up/down genes using z-score
get_all_outliers <- function(myMergeDF, getTop, thresh, comparison, cancer_genes) {
  
  # Filter in Patient: TPM > 10
  myMergeDF <- myMergeDF[myMergeDF$SampleX > 10,] 
  
  # z-score
  output <- apply(myMergeDF, FUN = get_zscore, MARGIN = 1) 
  
  # full data
  # combine tpm and z-score for sample of interest
  genes.df <- data.frame(z_score = output, tpm = myMergeDF[names(output),"SampleX"])
  genes.df$diff_expr <- ""
  genes.df$diff_expr[which(genes.df$z_score < (-1*thresh))] <- "down"
  genes.df$diff_expr[which(genes.df$z_score > thresh)] <- "up"
  genes.df$comparison <- comparison
  genes.df$cancer_gene <- ifelse(rownames(genes.df) %in% cancer_genes, TRUE, FALSE)
  
  # top 20 up/down genes 
  outputCanc <- output[intersect(names(output), cancer_genes)] # filter to cancer gene list
  outputUp <- sort(outputCanc, decreasing = T)[1:getTop] # top 20 up
  outputDown <- sort(outputCanc, decreasing = F)[1:getTop] # top 20 down
  outputUpDF <- data.frame(z_score = outputUp, tpm = myMergeDF[names(outputUp),"SampleX"])
  outputDownDF <- data.frame(z_score = outputDown, tpm = myMergeDF[names(outputDown),"SampleX"])
  diffexpr.top20 = rbind(outputUpDF, outputDownDF)
  
  return(list(expr.genes.z.score = output,
              diffexpr.top20 = diffexpr.top20,
              genes.df = genes.df))
}

# fgsea enrichment
fgsea_run <- function(pathways = gene_set, ranks = expr.genes.z.score, comparison){
  set.seed(42)
  fgseaRes <- fgsea(pathways, ranks, minSize = 15, maxSize = 1500, eps = 0.0)
  
  # significant pathways
  fgseaRes <- fgseaRes %>%
    filter(padj < 0.05) %>%
    mutate(direction = ifelse(ES > 0, "up", "down")) %>%
    arrange(padj) %>%
    mutate(genes = sapply(leadingEdge, paste, collapse=",")) %>%
    mutate(comparison = comparison)  %>%
    dplyr::select(-c(log2err, leadingEdge))
  
  return(fgseaRes)
}

# function to tabulate DE and Pathway results
run_rnaseq_analysis <- function(exp.data, refData, thresh = 2, gene_set = gene_set, comparison, single_sample = FALSE, sample_name = "SampleX") {
  
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
  geneAnalysisOut <- get_all_outliers(myMergeDF = mergeDF, getTop = 20, thresh = thresh, comparison = comparison, cancer_genes = cancer_genes$Gene_Symbol)
  
  # pathway enrichment: sort genes by z-score
  expr.genes.z.score <- geneAnalysisOut$expr.genes.z.score
  expr.genes.z.score <- expr.genes.z.score[order(expr.genes.z.score, decreasing = TRUE)]
  
  # significant pathways (adj. pvalue < 0.05)
  pathway.df <- fgsea_run(pathways = gene_set, ranks = expr.genes.z.score, comparison)

  # map of genes + pathways
  pathway.df.exp <- pathway.df %>%
    separate_rows(genes, convert = TRUE) %>%
    dplyr::select(pathway, genes)
  
  # create empty df just so you can map to gene.df below
  if(nrow(pathway.df.exp) == 0){
    pathway.df.exp <- data.frame(pathway = "", genes = "")
  }
  
  # map pathways to individual genes
  genes.df <- geneAnalysisOut$genes.df  %>%
    filter(diff_expr != "") %>%
    rownames_to_column("genes") %>%
    left_join(pathway.df.exp, by = "genes") %>%
    group_by(genes) %>%
    mutate(pathway = toString(pathway)) %>%
    unique() 
  
  if(single_sample == TRUE){
    colnames(mergeDF)[ncol(mergeDF)] <- sample_name
    finalOut <- list(pathways = pathway.df, 
                     genes = genes.df,
                     tpm = mergeDF[ncol(mergeDF)],
                     expr.genes.z.score = geneAnalysisOut$expr.genes.z.score,
                     diffexpr.top20 = geneAnalysisOut$diffexpr.top20)
  } else {
    finalOut <- list(pathways = pathway.df, 
                     genes = genes.df)
  }
  
  return(finalOut)
}
