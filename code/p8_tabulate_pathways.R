# Script to read in excel files from previous pnoc reports
# Capture upregulated pathways for all genomically similar patients

options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m"))
gc()

# function to merge excel sheets
merge.excel <- function(nm) {
  fname <- as.character(nm[[1]])
  x <- read.xlsx(file = fname, sheetIndex = 1)  # up pathways
  y <- read.xlsx(file = fname, sheetIndex = 3)  # down pathways
  x <- rbind(x, y)
  
  # highly significant up/down pathways only
  x <- x %>%
    filter(ADJ_P_VAL < 0.01) %>%
    dplyr::select(-c(Freq))
  if(nrow(x) > 1){
    x <- as.data.frame(x)
    return(x)
  }
}

# function to get upregulated pathways from top genomically similar patients  
tabulate_pathways <- function(allCor, numNeighbors = 20) {
  
  fname1 <- paste0(topDir, '/Summary/sig_pathways_gen_similar.txt')
  fname2 <- paste0(topDir, '/Summary/cnv_pathways.txt')
  if(!file.exists(fname1) | !file.exists(fname2)){
    # top 20 genomically similar PNOC008 patients 
    colnames(allCor)[1] <- "Correlation"
    patSamples <- allCor[1:numNeighbors,]
    patSamples <- patSamples[grep("PNOC", patSamples$sample_barcode),'sample_barcode']
    patSamples <- c(patSamples, sampleInfo$subjectID)
    
    # read all previous reports
    patPath <- list.files(path = 'data', pattern = "*.xlsx", recursive = TRUE, full.names = T)
    patPath <- grep('PNOC008-05-CHOP', patPath, invert = T, value = T) # remove CHOP report to avoid confusion
    patPath <- patPath[grep('PNOC008-',  patPath)]
    patPath <- data.frame(fname = patPath)
    patPath$sample_name <- gsub('.*/|_summary.xlsx', '', patPath$fname)
    patPath$sample_name <- gsub('-CHOP', '', patPath$sample_name)
    patPath$to_match <- gsub('-[0]+', '-', patPath$sample_name)
    
    # subset to genomically similar PNOC008 patients
    patPath <- patPath[which(patPath$to_match %in% patSamples),1:2]
    
    # now read up/down pathways sheet for each of these patients
    patPath <- plyr::ddply(.data = patPath, .variables = 'sample_name', .fun = function(x) merge.excel(x))
    
    # top 20 genomically similar PBTA
    pbtaSamples <- allCor[1:numNeighbors,]
    pbtaSamples <- pbtaSamples[grep("^BS_", pbtaSamples$sample_barcode),'sample_barcode']
    
    # read precalculated enrichment for PBTA vs GTEx, PBTA vs PBTA (HGG) and PBTA vs PBTA (All)
    gtexBrain <- readRDS('data/Reference/GSEA/PBTA_vs_GTExBrain.RDS')
    gtexBrain <- gtexBrain[pbtaSamples]
    gtexBrain <- plyr::ldply(gtexBrain, .fun = function(x) return(x[[1]]), .id = 'sample_name')
    pbtaAll <- readRDS('data/Reference/GSEA/PBTA_vs_PBTA.RDS')
    pbtaAll <- pbtaAll[pbtaSamples]
    pbtaAll <- plyr::ldply(pbtaAll, .fun = function(x) return(x[[1]]), .id = 'sample_name')
    pbtaHGG <- readRDS('data/Reference/GSEA/PBTA_vs_PBTAHGG.RDS')
    pbtaHGG <- pbtaHGG[pbtaSamples]
    pbtaHGG <- plyr::ldply(pbtaHGG, .fun = function(x) return(x[[1]]), .id = 'sample_name')
    pbtaPath <- rbind(pbtaAll, pbtaHGG, gtexBrain)

    # now merge patPath (PNOC008) and pbtaPath (PBTA)
    totalPath <- rbind(patPath, pbtaPath)
    
    # highly significant up/down pathways only (adj. pvalue < 0.01)
    # pathways which are seen misregulated in at least 80% of genomically similar samples
    totalPath <- totalPath %>%
      filter(ADJ_P_VAL < 0.01) %>% 
      dplyr::select(-c(SET_SIZE, NUM_GENES_INPUT, OVERLAP)) %>%
      mutate(P_VAL = scientific(P_VAL, digits = 3),
             ADJ_P_VAL = scientific(ADJ_P_VAL, digits = 3)) %>%
      group_by(Pathway, Comparison, Direction) %>%
      mutate(Sample.count.per.pathway = n()) %>%
      filter(Sample.count.per.pathway >= 16) %>%
      arrange(desc(sample_name), desc(Sample.count.per.pathway)) 
    
    # now create table2 in which we will have genes, pathway and copy number info
    # this is only for PNOC008 patient of interest
    # cnv gain/loss with WilcoxonRankSumTestPvalue < 0.05
    cnvPath <- cnvGenes %>%
      filter(Status %in% c("gain", "loss") & Pvalue < 0.05)
    cnvPath$Status <- ifelse(cnvPath$Status == "gain", "Amplification", "Deletion")
    cnvPath <- totalPath %>%
      filter(sample_name == sampleInfo$subjectID) %>%
      ungroup() %>%
      mutate(Pathway = paste0(Pathway,' (',Direction,')')) %>%
      dplyr::select(sample_name, Pathway, GENES, Direction, Comparison) %>%
      separate_rows(GENES) %>%
      filter(GENES != 1) %>%
      group_by(sample_name, GENES, Comparison) %>%
      summarise(Pathway = toString(Pathway)) %>%
      inner_join(cnvPath, by = c("GENES"="Gene"))
    
    # write out both tables
    write.table(x = totalPath, file = fname1, quote = F, sep = "\t", row.names = F)
    write.table(x = cnvPath, file = fname2, quote = F, sep = "\t", row.names = F)
  } else {
    totalPath <- read.delim(fname1, check.names = F)
    cnvPath <- read.delim(fname2, check.names = F)
  }
  return(list(totalPath, cnvPath))
}