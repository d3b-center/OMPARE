# Author: Komal S. Rathi
# Function: Connectivity Analysis using LINCS

suppressPackageStartupMessages(library(signatureSearch))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ExperimentHub))
suppressPackageStartupMessages(library(rhdf5))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(HDF5Array))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(visNetwork))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(ggpubr))

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))
dsigdb_dir <- file.path(ref_dir, "dsigdb")

# get LINCS data
eh <- ExperimentHub()
lincs <- eh[["EH3226"]]

# function to save network to file
network_to_file <- function(dnet_object, filename){
  ## Extract nodes from dtnetplot list object and remove shape parameter
  vertices <- dnet_object$x$nodes
  vertices$shape <- NULL
  
  ## Convert edge and node data frames to graph object
  net <- graph_from_data_frame(d = dnet_object$x$edges, vertices = vertices, directed = F)
  
  ## Assign layout
  l <- layout_with_fr(net)
  l <- norm_coords(l, ymin = -1, ymax = 1, xmin = -1, xmax = 1)
  pdf(file = filename, width = 10, height = 10)
  plot(net, layout = l*2, vertex.size = 4, rescale = T, vertex.label.cex = 0.5)
  dev.off()
}


drug_barplots <- function(dat, xlab, ylab, top = 20, fill_var = NULL, title){
  dat <- dat %>%
    dplyr::select(xlab, ylab, fill_var) %>%
    unique() %>%
    filter(get(ylab) != 0) %>%
    arrange(get(ylab)) %>%
    slice_head(n = top) %>%
    as.data.frame()
  
  dat[,xlab] <- factor(dat[,xlab], levels = unique(dat[,xlab]))
  if(!is.null(fill_var)){
    p <- ggplot(dat, aes(reorder(get(xlab), -get(ylab)), y = (-1)*log10(get(ylab)), fill = get(fill_var))) + 
      geom_bar(stat="identity") + coord_flip() + theme_bw() +
      xlab("") + 
      ylab("-log10 Adj. P-Value") + 
      scale_fill_manual(name = "Direction", values = c("down" = "forest green", "up" = "red")) +
      theme(plot.margin = unit(c(1, 5, 1, 7), "cm")) + 
      ggtitle(title)
  } else {
    p <- ggplot(dat, aes(get(xlab), y = (-1)*log10(get(ylab)))) + 
      geom_bar(stat="identity") + coord_flip() + theme_bw() +
      xlab("") + 
      ylab("-log10 Adj. P-Value") + 
      theme(plot.margin = unit(c(1, 5, 1, 7), "cm")) + 
      ggtitle(title) +
      scale_x_discrete(labels = function(x) str_wrap(x, width = 50))
  }
  
  return(p)
}

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
    drugs <- qSig_output %>%
      filter(trend == trend_val & WTCS_FDR < wtcs_fdr_cutoff) %>%
      .$pert %>%
      unique()
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
    drugs <- qSig_output %>%
      filter(cor_score < cor_score_cutoff) %>%
      .$pert %>%
      unique()
  }
  fname <- file.path(output_dir, paste0(comparison, "_qSig_output.txt"))
  write.table(qSig_output, file = fname, quote = F, sep = "\t", row.names = F) 
  
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
                      top = 20, fill_var = "trend",
                      title = "Query Signature")
  p2 <- drug_barplots(dat = tsea_reactome_df, 
                      xlab = "Description", ylab = "p.adjust",
                      top = 20, fill_var = NULL,
                      title = "TSEA Reactome")
  p3 <- drug_barplots(dat = dsea_go_mf_df, 
                      xlab = "Description", ylab = "p.adjust",
                      top = 20, fill_var = NULL,
                      title = "DSEA GO MF")
  p <- plot_grid(p1, p2, p3, ncol = 1, nrow = 3, axis = "lr", align = "v")
  return(p)
}

