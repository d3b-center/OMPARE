# Author: Komal S. Rathi
# Function: Connectivity Analysis using LINCS

suppressPackageStartupMessages({
  library(signatureSearch)
  library(tidyverse)
  library(dplyr)
  library(optparse)
  library(ExperimentHub)
  library(rhdf5)
  library(SummarizedExperiment)
  library(HDF5Array)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(visNetwork)
  library(igraph)
  library(cowplot)
  library(gridExtra)
  library(ggpubr)
})

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")

# source functions
source(file.path(module_dir, "utils", "network_to_file.R"))
source(file.path(module_dir, "utils", "drug_barplots.R"))

# get LINCS data
eh <- ExperimentHub()
lincs <- eh[["EH3226"]]

# MOA terms to drug name mappings obtained from Touchstone database at CLUE website 
data('clue_moa_list')
touchstone_data <- unlist(clue_moa_list)

# run connectivity analysis on each cluster
lincs_connectivity <- function(input, num_features = 2000, num_sets = 25, method = c("LINCS", "Cor"), wtcs_fdr_cutoff = 0.05, trend_val = c("up", "down"), cor_score_cutoff = 0, output_dir){
  
  comparison <- unique(gsub('_[0-9].*', '', input$comparison))
  
  # up and down genes
  cluster_upset <- input %>%
    filter(diff_expr == "up") %>%
    arrange(logFC) %>%
    slice_head(n = num_features)
  cluster_downset <- input %>%
    filter(diff_expr == "down") %>%
    arrange(logFC) %>%
    slice_head(n = num_features)
  
  if(method == "LINCS"){
    # map to ENTREZ identifiers
    upset = mapIds(org.Hs.eg.db, keys = cluster_upset$genes, column = "ENTREZID", keytype = "SYMBOL")
    downset = mapIds(org.Hs.eg.db, keys = cluster_downset$genes, column = "ENTREZID", keytype = "SYMBOL")
    
    # LINCS-based similarity metric
    qSig_output <- qSig(query = list(upset = upset, downset = downset), gess_method="LINCS", refdb = lincs)
    qSig_output <- gess_lincs(qSig = qSig_output, sortby = "NCS", tau = T, workers = 1)
    qSig_output <- result(qSig_output)
    
    # filter drugs
    qSig_output <- qSig_output %>%
      filter(trend == trend_val & WTCS_FDR < wtcs_fdr_cutoff,
             pert %in% touchstone_data)
    drugs <- unique(qSig_output$pert)
  } else if(method == "Cor"){
    # correlation-based similarity
    query_mat <- cluster_upset %>%
      rbind(cluster_downset) %>% 
      mutate(id = mapIds(org.Hs.eg.db, keys = geneSymbol, column = "ENTREZID", keytype = "SYMBOL")) %>% 
      arrange(id) %>% 
      filter(!is.na(id)) %>%
      column_to_rownames('id') %>%
      dplyr::select(score) %>% as.matrix()
    qSig_output <- qSig(query = query_mat, gess_method = "Cor", refdb = lincs)
    qSig_output <- gess_cor(qSig = qSig_output, method = "spearman", workers = 1)
    qSig_output <- result(qSig_output)
    
    # filter drugs
    qSig_output <- qSig_output %>%
      filter(cor_score < cor_score_cutoff,
             pert %in% touchstone_data)
    drugs <- unique(qSig_output$pert)
  }
  fname <- file.path(output_dir, paste0(comparison, "_qSig_output.txt"))
  write.table(qSig_output, file = fname, quote = F, sep = "\t", row.names = F) 
  
  if(length(drugs) != 0){
    # hypergeometric TSEA using Reactome 
    tsea_reactome <- tsea_dup_hyperG(drugs = drugs, 
                                     type = "Reactome", 
                                     pvalueCutoff = 0.5, 
                                     dt_anno = 'DrugBank',
                                     qvalueCutoff = 0.5, readable = TRUE)
    tsea_reactome_df <- result(tsea_reactome)
    fname <- file.path(output_dir, paste0(comparison, "_tsea_reactome_output.txt"))
    write.table(tsea_reactome_df, file = fname, quote = F, sep = "\t", row.names = F) 
    
    # hypergeometric DSEA using GO MF
    dsea_go_mf <- dsea_hyperG(drugs = drugs, type = "GO", ont = "MF")
    dsea_go_mf_df <- result(dsea_go_mf)
    fname <- file.path(output_dir, paste0(comparison, "_dsea_go_mf_output.txt"))
    write.table(dsea_go_mf_df, file = fname, quote = F, sep = "\t", row.names = F)
    
    # network plot
    top_set <- dsea_go_mf_df %>% arrange(p.adjust) %>% slice_head(n = num_sets) %>% .$ID
    dnet_object <- dtnetplot(drugs = drugs(dsea_go_mf), set = top_set, ont = "MF")
    fname <- file.path(output_dir, paste0(comparison, "_dsea_go_mf_output.html"))
    visSave(dnet_object, fname)
    fname <- file.path(output_dir, paste0(comparison, "_dsea_go_mf_output.pdf"))
    network_to_file(dnet_object = dnet_object, filename = fname)
    
    # plots
    p1 <- drug_barplots(dat = qSig_output, 
                        xlab = "pert", ylab = "WTCS_FDR",
                        top = 20, fill_var = NULL,
                        title = "Query Signature")
    p2 <- drug_barplots(dat = tsea_reactome_df, 
                        xlab = "Description", ylab = "p.adjust",
                        top = 20, fill_var = NULL,
                        title = "TSEA Reactome")
    p3 <- drug_barplots(dat = dsea_go_mf_df, 
                        xlab = "Description", ylab = "p.adjust",
                        top = 20, fill_var = NULL,
                        title = "DSEA GO MF")
    title <- cowplot::ggdraw() + draw_label(comparison, fontface = 'bold')
    p <- cowplot::plot_grid(title, 
              plot_grid(p1, p2, p3, ncol = 3, axis = "lr", align = "h"), 
              ncol = 1, rel_heights = c(0.05, 1))
  } else {
    p <- ggplot()
  }

  return(p)
}

