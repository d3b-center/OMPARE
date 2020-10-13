# capture upregulated pathways for all genomically similar patients

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))
gsea.dir <- file.path(ref_dir, 'gsea')

# function to get upregulated pathways from top genomically similar patients  
pathway_analysis <- function(all_cor, prefix, comparison) {
  
  # top 20 genomically similar PNOC008 patients 
  pnoc008_samples <- all_cor[grep("PNOC", all_cor$nearest_neighbor),'nearest_neighbor']
  pnoc008_samples <- c(pnoc008_samples, sampleInfo$subjectID)
  
  # 008 comparisons
  pnoc008_vs_gtex_brain <- readRDS(file.path(gsea.dir, 'pnoc008_vs_gtex_brain.rds'))
  pnoc008_vs_gtex_brain <- pnoc008_vs_gtex_brain[pnoc008_samples]
  pnoc008_vs_gtex_brain <- plyr::ldply(pnoc008_vs_gtex_brain, .fun = function(x) return(x[[1]]), .id = 'sample_name')
  pnoc008_vs_pbta <- readRDS(file.path(gsea.dir, 'pnoc008_vs_pbta.rds'))
  pnoc008_vs_pbta <- pnoc008_vs_pbta[pnoc008_samples]
  pnoc008_vs_pbta <- plyr::ldply(pnoc008_vs_pbta, .fun = function(x) return(x[[1]]), .id = 'sample_name')
  pnoc008_vs_pbta_hgg <- readRDS(file.path(gsea.dir, 'pnoc008_vs_pbta_hgg.rds'))
  pnoc008_vs_pbta_hgg <- pnoc008_vs_pbta_hgg[pnoc008_samples]
  pnoc008_vs_pbta_hgg <- plyr::ldply(pnoc008_vs_pbta_hgg, .fun = function(x) return(x[[1]]), .id = 'sample_name')
  pnoc008_vs_tcga_gbm <- readRDS(file.path(gsea.dir, 'pnoc008_vs_tcga_gbm.rds'))
  pnoc008_vs_tcga_gbm <- pnoc008_vs_tcga_gbm[pnoc008_samples]
  pnoc008_vs_tcga_gbm <- plyr::ldply(pnoc008_vs_tcga_gbm, .fun = function(x) return(x[[1]]), .id = 'sample_name')
  pnoc008_comparisons <- rbind(pnoc008_vs_pbta, pnoc008_vs_pbta_hgg, pnoc008_vs_gtex_brain, pnoc008_vs_tcga_gbm)
  
  # read precalculated enrichment for pbta vs gtex brain, pbta vs pbta hgg and pbta vs pbta
  if(comparison == "pediatric"){
    # top 20 genomically similar pbta
    pbta_samples <- all_cor[grep("^BS_", all_cor$nearest_neighbor),'nearest_neighbor']

    # comparisons
    pbta_vs_gtex_brain <- readRDS(file.path(gsea.dir, 'pbta_vs_gtex_brain.rds'))
    pbta_vs_gtex_brain <- pbta_vs_gtex_brain[pbta_samples]
    pbta_vs_gtex_brain <- plyr::ldply(pbta_vs_gtex_brain, .fun = function(x) return(x[[1]]), .id = 'sample_name')
    pbta_vs_pbta <- readRDS(file.path(gsea.dir, 'pbta_vs_pbta.rds'))
    pbta_vs_pbta <- pbta_vs_pbta[pbta_samples]
    pbta_vs_pbta <- plyr::ldply(pbta_vs_pbta, .fun = function(x) return(x[[1]]), .id = 'sample_name')
    pbta_vs_pbta_hgg <- readRDS(file.path(gsea.dir, 'pbta_vs_pbta_hgg.rds'))
    pbta_vs_pbta_hgg <- pbta_vs_pbta_hgg[pbta_samples]
    pbta_vs_pbta_hgg <- plyr::ldply(pbta_vs_pbta_hgg, .fun = function(x) return(x[[1]]), .id = 'sample_name')
    other_comparisons <- rbind(pbta_vs_pbta, pbta_vs_pbta_hgg, pbta_vs_gtex_brain)
  } else {
    # top 20 genomically similar tcga
    tcga_samples <- all_cor[grep("^TCGA", all_cor$nearest_neighbor),'nearest_neighbor']

    # comparisons
    tcga_gbm_vs_gtex_brain <- readRDS(file.path(gsea.dir, 'tcga_gbm_vs_gtex_brain.rds'))
    tcga_gbm_vs_gtex_brain <- tcga_gbm_vs_gtex_brain[tcga_samples]
    tcga_gbm_vs_gtex_brain <- plyr::ldply(tcga_gbm_vs_gtex_brain, .fun = function(x) return(x[[1]]), .id = 'sample_name')
    tcga_gbm_vs_tcga_gbm <- readRDS(file.path(gsea.dir, 'tcga_gbm_vs_tcga_gbm.rds'))
    tcga_gbm_vs_tcga_gbm <- tcga_gbm_vs_tcga_gbm[tcga_samples]
    tcga_gbm_vs_tcga_gbm <- plyr::ldply(tcga_gbm_vs_tcga_gbm, .fun = function(x) return(x[[1]]), .id = 'sample_name')
    other_comparisons <- rbind(tcga_gbm_vs_gtex_brain, tcga_gbm_vs_tcga_gbm)
  }
  
  # now combine pnoc008_comparisons and other_comparisons
  shared_pathways <- rbind(pnoc008_comparisons, other_comparisons)
  
  # different cutoff for pediatric and adult tumors
  if(comparison == "pediatric"){
    n.perc = 18 # 90%
  } else {
    n.perc = 12 # 60%
  }

  # highly significant up/down pathways only (adj. pvalue < 0.01)
  # pathways which are seen misregulated in n.perc genomically similar samples
  shared_pathways <- shared_pathways %>%
    filter(padj < 0.01) %>% 
    dplyr::select(-c(ES, NES, size)) %>%
    mutate(pval = scientific(pval, digits = 3),
           padj = scientific(padj, digits = 3)) %>%
    group_by(pathway, comparison, direction) %>%
    mutate(Sample.count.per.pathway = n()) %>%
    filter(Sample.count.per.pathway >= n.perc) %>%
    arrange(desc(sample_name), desc(Sample.count.per.pathway)) %>%
    as.data.frame()

  # now create table2 in which we will have genes, pathway and copy number info
  # this is only for PNOC008 patient of interest
  # cnv gain/loss with WilcoxonRankSumTestPvalue < 0.05
  if(comparison == "pediatric"){
    cnv_mapping <- cnvGenes %>%
      filter(status %in% c("gain", "loss") & WilcoxonRankSumTestPvalue < 0.05)
    cnv_mapping$status <- ifelse(cnv_mapping$status == "gain", "Amplification", "Deletion")
    cnv_mapping <- shared_pathways %>%
      filter(sample_name == sampleInfo$subjectID) %>%
      ungroup() %>%
      mutate(pathway = paste0(pathway,' (', direction,')')) %>%
      dplyr::select(sample_name, pathway, genes, direction, comparison) %>%
      separate_rows(genes) %>%
      filter(genes != 1) %>%
      group_by(sample_name, genes, comparison) %>%
      summarise(pathway = toString(pathway)) %>%
      inner_join(cnv_mapping, by = c("genes"="hgnc_symbol"))
  } else {
    cnv_mapping <- data.frame()
  }
  return(list(shared_pathways = shared_pathways, cnv_mapping = cnv_mapping))
}
