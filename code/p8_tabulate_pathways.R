# Script to read in excel files from previous pnoc reports
# Capture upregulated pathways for all genomically similar patients

options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m"))

# function to merge excel sheets
merge.excel <- function(nm) {
  fname <- as.character(nm[[1]])
  x <- read.xlsx(file = fname, sheetIndex = 1)
  x <- x %>%
    mutate(Pathway.count.per.sample = Freq,
           P_VAL = scientific(P_VAL, digits = 3),
           ADJ_P_VAL = scientific(ADJ_P_VAL, digits = 3)) %>%
    dplyr::select(Pathway, Comparison, Direction, P_VAL, ADJ_P_VAL, Pathway.count.per.sample)
  if(nrow(x) > 1){
    x <- as.data.frame(x)
    return(x)
  }
}

# function to get upregulated pathways from top genomically similar patients  
tabulate_pathways <- function(allCor, numNeighbors = 20) {
  
  fname <- paste0(topDir, '/Summary/up_pathways_gen_similar.txt')
  if(!file.exists(fname)){
    # top 20 genomically similar pnoc patients 
    colnames(allCor)[1] <- "Correlation"
    patSamples <- allCor[1:numNeighbors,]
    patSamples <- patSamples[grep("PNOC", patSamples$sample_barcode),'sample_barcode']
    patSamples <- c(patSamples, sampleInfo$subjectID)
    
    # all previous reports
    pat.expDat <- list.files(path = 'data', pattern = "*.xlsx", recursive = TRUE, full.names = T)
    pat.expDat <- grep('PNOC008-05-NANT', pat.expDat, invert = T, value = T) # remove NANT report to avoid confusion
    pat.expDat <- pat.expDat[grep('PNOC008-',  pat.expDat)]
    pat.expDat <- data.frame(fname = pat.expDat)
    pat.expDat$sample_name <- gsub('.*/|_summary.xlsx', '', pat.expDat$fname)
    pat.expDat$sample_name <- gsub('-CHOP', '', pat.expDat$sample_name)
    pat.expDat$to_match <- gsub('-[0]+', '-', pat.expDat$sample_name)
    
    # subset to genomically similar patients
    pat.expDat <- pat.expDat[which(pat.expDat$to_match %in% patSamples),1:2]
    
    # now read upregulated pathways sheet for each of these patients
    pat.expr.mat <- plyr::ddply(.data = pat.expDat, .variables = 'sample_name', .fun = function(x) merge.excel(x))
    pat.expr.mat <- pat.expr.mat %>%
      group_by(Pathway, Comparison) %>%
      mutate(Sample.count.per.pathway = n()) %>%
      arrange(desc(sample_name), desc(Sample.count.per.pathway))
    write.table(x = pat.expr.mat, file = fname, quote = F, sep = "\t", row.names = F)
  } else {
    pat.expr.mat <- read.delim(fname, check.names = F)
  }
  return(pat.expr.mat)
}
