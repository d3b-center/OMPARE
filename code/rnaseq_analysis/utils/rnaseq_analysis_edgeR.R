# Author: Komal S. Rathi
# fgsea based pathway enrichment

suppressPackageStartupMessages({
  library(tidyverse)
  library(plyr)
  library(fgsea)
})

# source functions for edgeR's single sample differential analysis
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "rnaseq_analysis")

# source normalization function
source(file.path(module_dir, "utils", "get_all_outliers.R"))
source(file.path(module_dir, "utils", "ss_diffexpr.R"))
source(file.path(module_dir, "utils", "fgsea_run.R"))

# function to tabulate DE and Pathway results
run_rnaseq_analysis_edger <- function(exp.data.counts, exp.data.tpm, refData.counts, gene_set = gene_set, comparison, single_sample = FALSE, sample_name = "SampleX", cancer_genes) {
  
  # add genes as rownames and remove sample name
  exp.data.counts <- exp.data.counts %>%
    remove_rownames() %>%
    column_to_rownames('gene_symbol') %>%
    dplyr::select(-c(sample))
  
  # now remove the current sample from refData (required when we are doing sample from the same cohort)
  refData.counts <- refData.counts[,grep(sample_name, colnames(refData.counts), invert = T, value = T)]
  
  # Merge refData and sample of interest data on common genes
  mergeDF.counts <- refData.counts %>% rownames_to_column('gene_symbol') %>% 
    inner_join(exp.data.counts %>% 
                 rownames_to_column('gene_symbol'), by = 'gene_symbol') %>% 
    column_to_rownames('gene_symbol')
  colnames(mergeDF.counts)[ncol(mergeDF.counts)] <- "sample_of_interest" # sample of interest
  
  # subset tpm data with same gene symbols
  exp.data.tpm <- exp.data.tpm  %>%
    filter(sample %in% sample_name) %>%
    filter(gene_symbol %in% rownames(mergeDF.counts))
  
  # get gene outliers (top 20 Up and Down)
  geneAnalysisOut <- get_all_outliers(myMergeDF.counts = mergeDF.counts, 
                                      expr.tpm = exp.data.tpm,
                                      getTop = 20, 
                                      comparison = comparison, 
                                      cancer_genes = cancer_genes$Gene_Symbol)
  
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
      dplyr::mutate(!!sample_name := tpm) %>%
      dplyr::select(-c(sample, tpm)) %>%
      column_to_rownames('gene_symbol')
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
