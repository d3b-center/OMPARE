# Author: Komal S. Rathi
# fgsea based pathway enrichment

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(fgsea))

# source functions for edgeR's single sample differential analysis
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(patient_level_analyses_utils, 'rnaseq_edger_normalizations.R'))

# housekeeping genes (optional)
housekeeping_genes <- read.csv(file.path(ref_dir, 'Housekeeping_GenesHuman.csv'), sep = ';')
housekeeping_genes <- unique(housekeeping_genes$Gene.name)

# get up/down genes using edgeR
get_all_outliers <- function(myMergeDF.counts, expr.tpm, getTop, thresh, comparison, cancerGeneNames){
  
  genes.df <- ss_expr(expr = myMergeDF.counts, norm_method = "tmm", housekeeping_genes = NULL)
  genes.df <- genes.df %>%
    rownames_to_column('Gene') %>%
    inner_join(expr.tpm, by = 'Gene') %>%
    mutate(tpm = TPM) %>%
    dplyr::select(logFC, tpm, diff_expr, Gene) %>%
    column_to_rownames('Gene')
  
  # add comparison and annotate cancer genes 
  genes.df$comparison <- comparison
  genes.df$cancer_gene <- ifelse(rownames(genes.df) %in% cancerGeneNames, TRUE, FALSE)
  
  # top 20 up/down genes 
  output <- genes.df$logFC
  names(output) <- rownames(genes.df)
  outputCanc <- output[intersect(names(output), cancerGeneNames)] # filter to cancer gene list
  outputUp <- sort(outputCanc, decreasing = T)[1:getTop] # top 20 up
  outputDown <- sort(outputCanc, decreasing = F)[1:getTop] # top 20 down
  outputUpDF <- data.frame(logfc = outputUp, tpm = genes.df[names(outputUp),"tpm"])
  outputDownDF <- data.frame(logfc = outputDown, tpm = genes.df[names(outputDown),"tpm"])
  diffexpr.top20 = rbind(outputUpDF, outputDownDF)
  
  return(list(expr.genes.logfc = output,
              diffexpr.top20 = diffexpr.top20,
              genes.df = genes.df))
}

# fgsea enrichment
fgsea_run <- function(pathways, ranks, comparison){
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
run_rnaseq_analysis_edger <- function(exp.data.counts, exp.data.tpm, refData.counts, gene_set = gene_set, comparison, single_sample = FALSE, sample_name = "SampleX") {
  
  # current sample of interest
  smp <- unique(as.character(exp.data.counts$Sample))
  
  # add genes as rownames and remove sample name
  exp.data.counts <- exp.data.counts %>%
    remove_rownames() %>%
    column_to_rownames('Gene') %>%
    dplyr::select(-c(Sample))
  
  # now remove the current sample from refData (required when we are doing sample from the same cohort)
  refData.counts <- refData.counts[,grep(smp, colnames(refData.counts), invert = T, value = T)]
  
  # Merge refData and sample of interest data on common genes
  mergeDF.counts <- refData.counts %>% rownames_to_column('gene') %>% 
    inner_join(exp.data.counts %>% 
                 rownames_to_column('gene'), by = 'gene') %>% 
    column_to_rownames('gene')
  colnames(mergeDF.counts)[ncol(mergeDF.counts)] <- "SampleX" # sample of interest
  
  # subset tpm data with same gene symbols
  exp.data.tpm <- exp.data.tpm  %>%
    filter(Sample %in% smp) %>%
    filter(Gene %in% rownames(mergeDF.counts))
  
  # get gene outliers (top 20 Up and Down)
  geneAnalysisOut <- get_all_outliers(myMergeDF.counts = mergeDF.counts, 
                                      expr.tpm = exp.data.tpm,
                                      getTop = 20, 
                                      comparison = comparison, 
                                      cancerGeneNames = cancer_genes$Gene_Symbol)
  
  # pathway enrichment: sort genes by logFC
  expr.genes.logfc <- geneAnalysisOut$expr.genes.logfc
  expr.genes.logfc <- expr.genes.logfc[order(expr.genes.logfc, decreasing = TRUE)]
  
  # significant pathways (adj. pvalue < 0.05)
  pathway.df <- fgsea_run(pathways = gene_set, ranks = expr.genes.logfc, comparison)
  
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
    exp.data.tpm <- exp.data.tpm %>%
      mutate(!!sample_name := TPM) %>%
      dplyr::select(!!sample_name, Gene) %>%
      column_to_rownames('Gene')
    finalOut <- list(pathways = pathway.df, 
                     genes = genes.df,
                     tpm = exp.data.tpm,
                     expr.genes.logfc = geneAnalysisOut$expr.genes.logfc,
                     diffexpr.top20 = geneAnalysisOut$diffexpr.top20)
  } else {
    finalOut <- list(pathways = pathway.df, 
                     genes = genes.df)
  }
  
  return(finalOut)
}
