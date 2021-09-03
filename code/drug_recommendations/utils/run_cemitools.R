suppressPackageStartupMessages({
  library(CEMiTool)
  library(tidyverse)
  library(patchwork)
  library(dplyr)
  library(callr)
})

# function to run cemitools analysis
run_cemitools <- function(exp_counts_corrected, clinical, output_dir, cemitools_dir, pnoc008_cluster){
  # run CEMItools
  n <- 100
  cem <- cemitool(as.data.frame(exp_counts_corrected), 
                  clinical, 
                  filter = T,
                  cor_functio = 'bicor', 
                  network_type = 'signed',
                  tom_type = 'signed',
                  sample_name_column = 'Kids_First_Biospecimen_ID',
                  class_column = 'CC',
                  merge_similar = T,
                  apply_vst = T, 
                  verbose = F)
  
  # write out hubs and summary
  hubs <- get_hubs(cem, n, method = "kME")
  summary <- mod_summary(cem)
  saveRDS(hubs, file = file.path(cemitools_dir, 'hubs.rds'))
  saveRDS(summary, file = file.path(cemitools_dir, 'summary.rds'))
  
  # generate heatmap of gene set enrichment analysis
  cem <- mod_gsea(cem)
  cem <- plot_gsea(cem)
  
  # plot gene expression within each module
  cem <- plot_profile(cem)
  
  # read GMT file - reactome file
  gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
  gmt_in <- CEMiTool::read_gmt(gmt_fname)
  
  # perform over representation analysis
  cem <- mod_ora(cem, gmt_in)
  
  # plot ora results
  cem <- plot_ora(cem)
  
  # read interactions
  int_fname <- system.file("extdata", "interactions.tsv", package = "CEMiTool")
  int_df <- read.delim(int_fname)
  
  # plot interactions
  interactions_data(cem) <- int_df # add interactions
  cem <- plot_interactions(cem) # generate plot
  
  # save interaction plots for all modules
  r(function(x, y) { CEMiTool::diagnostic_report(cem = x, directory = y, force = T) }, 
    args = list(cem, cemitools_dir))

  # output
  r(function(x, y) { CEMiTool::generate_report(cem = x, directory = y, force = T) }, 
    args = list(cem, cemitools_dir))
  
  write_files(cem, directory = cemitools_dir, force = T)
  save_plots(cem, "all", directory = cemitools_dir, force = T)
  
  # get pos/neg correlated modules for pnoc008 cluster using p-adj < 0.05 cut-off
  corr_modules <- cem@enrichment$padj %>%
    filter(get(as.character(pnoc008_cluster)) < 0.05) %>%
    .$pathway
  corr_modules <- cem@enrichment$nes %>%
    filter(pathway %in% corr_modules,
           pathway != "Not.Correlated") %>%
    mutate(direction = ifelse(get(as.character(pnoc008_cluster)) > 0, "pos", "neg"))
  
  # get hub genes for pos/neg correlated modules with a cutoff of 0.5
  # note: when using method = "kME", the hubs object does not behave like a normal list and so I am unable to use stack to unlist the list recursively
  network_hubs <- stack(unlist(hubs, use.names = T))
  network_hubs$ind <- gsub('Not.Correlated','Not_Correlated', network_hubs$ind)
  network_hubs <-  cbind(values = network_hubs$values, reshape2::colsplit(network_hubs$ind, pattern = '\\.', names = c("module", "genes")))
  network_hubs <- network_hubs %>%
    filter(module %in% corr_modules$pathway,
           values >= 0.5)
  
  # annotate targetable hubs
  fname <- file.path(output_dir, "transcriptome_drug_rec.rds")
  dge_genes <- readRDS(fname)
  dge_genes <- dge_genes %>% 
    mutate(Network_Hub = ifelse(Comparison == "PBTA_HGG_189" & Gene %in% network_hubs$genes, "Yes", "No"))
  
  # rewrite transcriptiomic based drug recommendations
  saveRDS(dge_genes, file = fname)
  
  # get ora data for pos/neg correlated modules
  ora_dat <- ora_data(cem = cem)
  ora_dat <- ora_dat %>%
    inner_join(corr_modules %>% 
                 dplyr::select(pathway, direction), by = c("Module" = "pathway"))
  ora_dat <- ora_dat %>%
    group_by(Module) %>%
    dplyr::mutate(order_val = row_number()) %>%
    arrange(p.adjust) %>%
    slice_head(n = 10) 
  modules <- unique(ora_dat$Module)
  pdf(file.path(output_dir, "ora_plots.pdf"), width = 12, height = 8)
  for(i in 1:length(modules)){
    tmp <- ora_dat %>%
      filter(Module %in% modules[i])
    title <- paste("Module: ", unique(tmp$Module), "| Direction: ", unique(tmp$direction))
    p <- ggplot(tmp, aes(x = reorder(ID, -p.adjust), 
                         y = -log10(p.adjust), 
                         fill = -log10(p.adjust))) +
      geom_bar(stat = "identity") +
      coord_flip() +
      xlab('') + ylab('−log10(adjusted p−value)') + ggtitle(title) +
      theme_bw() +
      theme_Publication(base_size = 12) 
    print(p)
  }
  dev.off()
}