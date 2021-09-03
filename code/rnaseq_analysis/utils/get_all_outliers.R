# get up/down genes using edgeR
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "rnaseq_analysis")

# source function
source(file.path(module_dir, "utils", "ss_diffexpr.R"))

get_all_outliers <- function(myMergeDF.counts, expr.tpm, getTop, thresh, comparison, cancer_genes){
  
  genes.df <- ss_diffexpr(expr = myMergeDF.counts, norm_method = "tmm", housekeeping_genes = NULL)
  genes.df <- genes.df %>%
    rownames_to_column('gene_symbol') %>%
    inner_join(expr.tpm, by = 'gene_symbol') %>%
    dplyr::select(logFC, tpm, diff_expr, gene_symbol) %>%
    column_to_rownames('gene_symbol')
  
  # add comparison and annotate cancer genes 
  genes.df$comparison <- comparison
  genes.df$cancer_gene <- ifelse(rownames(genes.df) %in% cancer_genes, TRUE, FALSE)
  
  # top 20 up/down genes 
  output <- genes.df$logFC
  names(output) <- rownames(genes.df)
  outputCanc <- output[intersect(names(output), cancer_genes)] # filter to cancer gene list
  outputUp <- sort(outputCanc, decreasing = T)[1:getTop] # top 20 up
  outputDown <- sort(outputCanc, decreasing = F)[1:getTop] # top 20 down
  outputUpDF <- data.frame(logfc = outputUp, tpm = genes.df[names(outputUp),"tpm"])
  outputDownDF <- data.frame(logfc = outputDown, tpm = genes.df[names(outputDown),"tpm"])
  diffexpr.top20 = rbind(outputUpDF, outputDownDF)
  
  return(list(expr.genes.logfc = output,
              diffexpr.top20 = diffexpr.top20,
              genes.df = genes.df))
}