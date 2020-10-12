# capture upregulated pathways for all genomically similar patients

options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m"))
gc()

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))
gsea.dir <- file.path(ref_dir, 'GSEA')

# function to get upregulated pathways from top genomically similar patients  
pathway_analysis <- function(all_cor, prefix, comparison) {
  
  # top 20 genomically similar PNOC008 patients 
  patSamples <- all_cor[grep("PNOC", all_cor$nearest_neighbor),'nearest_neighbor']
  patSamples <- c(patSamples, sampleInfo$subjectID)
  
  # 008 comparisons
  gtexBrain <- readRDS(file.path(gsea.dir, 'PNOC008_vs_GTExBrain.RDS'))
  gtexBrain <- gtexBrain[patSamples]
  gtexBrain <- plyr::ldply(gtexBrain, .fun = function(x) return(x[[1]]), .id = 'sample_name')
  pbtaAll <- readRDS(file.path(gsea.dir, 'PNOC008_vs_PBTA.RDS'))
  pbtaAll <- pbtaAll[patSamples]
  pbtaAll <- plyr::ldply(pbtaAll, .fun = function(x) return(x[[1]]), .id = 'sample_name')
  pbtaHGG <- readRDS(file.path(gsea.dir, 'PNOC008_vs_PBTA_HGG.RDS'))
  pbtaHGG <- pbtaHGG[patSamples]
  pbtaHGG <- plyr::ldply(pbtaHGG, .fun = function(x) return(x[[1]]), .id = 'sample_name')
  tcgaGBM <- readRDS(file.path(gsea.dir, 'PNOC008_vs_TCGA_GBM.RDS'))
  tcgaGBM <- tcgaGBM[patSamples]
  tcgaGBM <- plyr::ldply(tcgaGBM, .fun = function(x) return(x[[1]]), .id = 'sample_name')
  patPath <- rbind(pbtaAll, pbtaHGG, gtexBrain, tcgaGBM)
  
  # read precalculated enrichment for PBTA vs GTEx, PBTA vs PBTA (HGG) and PBTA vs PBTA (All)
  if(comparison == "pediatric"){
    # top 20 genomically similar PBTA
    pbtaSamples <- all_cor[grep("^BS_", all_cor$nearest_neighbor),'nearest_neighbor']
    
    # comparisons
    gtexBrain <- readRDS(file.path(gsea.dir, 'PBTA_vs_GTExBrain.RDS'))
    gtexBrain <- gtexBrain[pbtaSamples]
    gtexBrain <- plyr::ldply(gtexBrain, .fun = function(x) return(x[[1]]), .id = 'sample_name')
    pbtaAll <- readRDS(file.path(gsea.dir, 'PBTA_vs_PBTA.RDS'))
    pbtaAll <- pbtaAll[pbtaSamples]
    pbtaAll <- plyr::ldply(pbtaAll, .fun = function(x) return(x[[1]]), .id = 'sample_name')
    pbtaHGG <- readRDS(file.path(gsea.dir, 'PBTA_vs_PBTAHGG.RDS'))
    pbtaHGG <- pbtaHGG[pbtaSamples]
    pbtaHGG <- plyr::ldply(pbtaHGG, .fun = function(x) return(x[[1]]), .id = 'sample_name')
    tumorPath <- rbind(pbtaAll, pbtaHGG, gtexBrain)
  } else {
    # top 20 genomically similar TCGA
    tcgaSamples <- all_cor[grep("^TCGA-", all_cor$nearest_neighbor),'nearest_neighbor']
    
    gtexBrain <- readRDS(file.path(gsea.dir, 'TCGA_GBM_vs_GTExBrain.RDS'))
    gtexBrain <- gtexBrain[tcgaSamples]
    gtexBrain <- plyr::ldply(gtexBrain, .fun = function(x) return(x[[1]]), .id = 'sample_name')
    tcgaGBM <- readRDS(file.path(gsea.dir, 'TCGA_GBM_vs_TCGA_GBM.RDS'))
    tcgaGBM <- tcgaGBM[tcgaSamples]
    tcgaGBM <- plyr::ldply(tcgaGBM, .fun = function(x) return(x[[1]]), .id = 'sample_name')
    tumorPath <- rbind(gtexBrain, tcgaGBM)
  }
  
  # now merge patPath (PNOC008) and tumorPath
  totalPath <- rbind(patPath, tumorPath)
  
  if(comparison == "pediatric"){
    n.perc = 18 # 90%
  } else {
    n.perc = 16 # 80%
  }
  # highly significant up/down pathways only (adj. pvalue < 0.01)
  # pathways which are seen misregulated in n.perc genomically similar samples
  totalPath <- totalPath %>%
    filter(ADJ_P_VAL < 0.01) %>% 
    dplyr::select(-c(SET_SIZE, NUM_GENES_INPUT, OVERLAP)) %>%
    mutate(P_VAL = scientific(P_VAL, digits = 3),
           ADJ_P_VAL = scientific(ADJ_P_VAL, digits = 3)) %>%
    group_by(Pathway, Comparison, Direction) %>%
    mutate(Sample.count.per.pathway = n()) %>%
    filter(Sample.count.per.pathway >= n.perc) %>%
    arrange(desc(sample_name), desc(Sample.count.per.pathway)) %>%
    as.data.frame()

  # now create table2 in which we will have genes, pathway and copy number info
  # this is only for PNOC008 patient of interest
  # cnv gain/loss with WilcoxonRankSumTestPvalue < 0.05
  if(comparison == "pediatric"){
    cnvPath <- cnvGenes %>%
      filter(status %in% c("gain", "loss") & WilcoxonRankSumTestPvalue < 0.05)
    cnvPath$status <- ifelse(cnvPath$status == "gain", "Amplification", "Deletion")
    cnvPath <- totalPath %>%
      filter(sample_name == sampleInfo$subjectID) %>%
      ungroup() %>%
      mutate(Pathway = paste0(Pathway,' (',Direction,')')) %>%
      dplyr::select(sample_name, Pathway, GENES, Direction, Comparison) %>%
      separate_rows(GENES) %>%
      filter(GENES != 1) %>%
      group_by(sample_name, GENES, Comparison) %>%
      summarise(Pathway = toString(Pathway)) %>%
      inner_join(cnvPath, by = c("GENES"="hgnc_symbol"))
  } else {
    cnvPath <- data.frame()
  }
  return(list(shared_pathways = totalPath, cnv_mapping = cnvPath))
}