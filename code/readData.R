#############################
# Read Patient-specific Data
# Save all data to global env
#############################

source('code/load_reference.R') # load reference data
source('code/filterFusions.R')  # load filter fusions script to filter fusions

readData <- function(topDir, fusion_method = c("star","arriba"), snv_pattern = "all"){
  
  # patient sample info (at most one with .txt extension)
  # assign n/a if no clinical info is present
  sampleInfo <- list.files(path = topDir, pattern = "patient_report.txt", recursive = T, full.names = T)
  if(length(sampleInfo) == 1){
    sampleInfo <- read.delim(sampleInfo, stringsAsFactors = F)
    assign("sampleInfo", sampleInfo, envir = globalenv())
  }
  
  # mutation data (can be multiple with .maf extension)
  # if snv_pattern is all, then use all files except consensus
  somatic.mut.pattern <- '*.maf'
  mutFiles <- list.files(path = topDir, pattern = somatic.mut.pattern, recursive = TRUE, full.names = T)
  if(snv_pattern == "all"){ 
    mutFiles <- grep('consensus', mutFiles, invert = TRUE, value = TRUE)
  } else {
    mutFiles <- grep(snv_pattern, mutFiles, value = TRUE)
  }
  if(length(mutFiles) >= 1){
    mutFiles <- lapply(mutFiles, data.table::fread, skip = 1, stringsAsFactors = F)
    mutData <- data.table::rbindlist(mutFiles)
    mutData <- as.data.frame(mutData)
    mutData <- unique(mutData)
    assign("mutData", mutData, envir = globalenv())
  } 
  
  # germline data
  mutFiles <- list.files(path = topDir, pattern = 'hg38_multianno.txt.gz', recursive = TRUE, full.names = T)
  if(length(mutFiles) >= 1){
    mutFiles <- lapply(mutFiles, data.table::fread, stringsAsFactors = F)
    mutData.germ <- data.table::rbindlist(mutFiles)
    mutData.germ <- as.data.frame(mutData.germ)
    mutData.germ <- unique(mutData.germ)
    assign("mutData.germ", mutData.germ, envir = globalenv())
  } 
  
  # copy number pvalue file
  cnvData <- list.files(path = topDir, pattern = "*.CNVs.p.value.txt", recursive = TRUE, full.names = T)
  if(length(cnvData) == 1){
    cnvData <- data.table::fread(cnvData, header = T, check.names = T)
    cnvData <- cnvData %>% 
      dplyr::select(chr, start, end, copy.number, 
                    status, WilcoxonRankSumTestPvalue) %>%
      as.data.frame()
    assign("cnvData", cnvData, envir = globalenv())
  }
  
  # copy number ratio data for plot
  cnvRatioData <- list.files(path = topDir, pattern = "*controlfreec.ratio.txt$", recursive = TRUE, full.names = T)
  if(length(cnvRatioData) == 1){
    cnvRatioData <- data.table::fread(cnvRatioData, header = T)
    cnvRatioData <- as.data.frame(cnvRatioData)
    assign("cnvRatioData", cnvRatioData, envir = globalenv())
  }
  
  # fusion data (chose either star or arriba or both)
  # fusion_method determines the pattern to be searched
  if(fusion_method == "star"){
    fusPattern = "*star-fusion.fusion_candidates.final"
  } else if(fusion_method == "arriba") {
    fusPattern = "*.arriba.fusions.tsv"
  } else {
    fusPattern = "*star-fusion.fusion_candidates.final|*.arriba.fusions.tsv"
  }
  # read fusion files + filter and merge them
  fusFiles <- list.files(path = topDir, pattern = fusPattern, recursive = TRUE, full.names = T)
  if(length(fusFiles) >= 1){
    fusFiles <- lapply(fusFiles, filterFusions)
    fusData <- do.call('rbind', fusFiles)
    fusData <- unique(fusData)
    if(nrow(fusData) >= 1){
      assign("fusData", fusData, envir = globalenv())
    }
  }
  
  # expression data: TPM (can be only 1 per patient with .genes.results)
  expDat <- list.files(path = topDir, pattern = "*.genes.results*", recursive = TRUE, full.names = T)
  if(length(expDat) == 1){
    expData <- read.delim(expDat)
    expData.full <- expData %>% 
      mutate(gene_id = str_replace(gene_id, "_PAR_Y_", "_"))  %>%
      separate(gene_id, c("gene_id", "gene_symbol"), sep = "\\_", extra = "merge") %>%
      mutate(gene_id = gsub('[.].*', '', gene_id))  %>%
      unique()
    expData <- expData.full %>% 
      arrange(desc(TPM)) %>% 
      distinct(gene_symbol, .keep_all = TRUE) %>%
      mutate(!!sampleInfo$subjectID := TPM) %>%
      dplyr::select(!!sampleInfo$subjectID, gene_id, gene_symbol) %>%
      unique()
    rownames(expData) <- expData$gene_symbol
    assign("expData", expData, envir = globalenv())
    
    expData.counts <- expData.full %>%
      filter(expData.full$gene_id  %in% expData$gene_id) %>%
      mutate(!!sampleInfo$subjectID := expected_count) %>%
      dplyr::select(!!sampleInfo$subjectID, gene_id, gene_symbol) %>%
      unique()
    rownames(expData.counts) <- expData.counts$gene_symbol
    assign("expData.counts", expData.counts, envir = globalenv())
  }
}

