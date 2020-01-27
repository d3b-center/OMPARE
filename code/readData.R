#############################
# Read Patient-specific Data
# Save all data to global env
#############################

source('code/helper.R')
source('code/load_reference.R')
source('code/filterFusions.R')

readData <- function(topDir, fusion_method = c("star","arriba"), snv_pattern = "all"){
  
  # patient sample info (at most one with .txt extension)
  # assign n/a if no clinical info is present
  sampleInfo <- list.files(path = paste0(topDir, "Clinical"), pattern = "*.txt", full.names = T)
  if(length(sampleInfo) == 1){
    sampleInfo <- read.delim(sampleInfo, stringsAsFactors = F)
  } else {
    sampleInfo <- setNames(data.frame(t(rep("n/a", 15)), stringsAsFactors = F), 
                           nm = c("patientName","reportDate","reportVersion","primRelapse",
                                  "tumorType","tumorLocation","collectionDate","labDirector",
                                  "pathologist","primPhysician","medicalFacility","ethnicity",
                                  "dob","sex","subjectID"))
  }
  assign("sampleInfo", sampleInfo, envir = globalenv())
  
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
  
  # also read mutect2 for all reports for TMB profile
  somatic.mut.pattern <- '*.maf'
  mutFiles <- list.files(path = topDir, pattern = somatic.mut.pattern, recursive = TRUE, full.names = T)
  mutFiles <- grep("mutect2", mutFiles, value = TRUE)
  if(length(mutFiles) >= 1){
    mutFiles <- lapply(mutFiles, data.table::fread, skip = 1, stringsAsFactors = F)
    mutect2Data <- data.table::rbindlist(mutFiles)
    mutect2Data <- as.data.frame(mutect2Data)
    mutect2Data <- unique(mutect2Data)
    assign("mutect2Data", mutect2Data, envir = globalenv())
  }
  
  # TMB scores
  TMBFiles <- list.files(path = topDir, pattern = "TMB", recursive = T, full.names = T)
  if(length(TMBFiles) == 2) {
    pedTMB <- grep('TCGA', TMBFiles, invert = T, value = T)
    pedTMB <- data.table::fread(pedTMB)
    assign("pedTMB", pedTMB, envir = globalenv())
    adultTMB <- grep('TCGA', TMBFiles, value = T)
    adultTMB <- data.table::fread(adultTMB)
    assign("adultTMB", adultTMB, envir = globalenv())
  }
  
  # TMB bedfile
  TMBFileBED <- list.files(path = topDir, pattern = ".bed$", recursive = T, full.names = T)
  if(length(TMBFileBED) == 1) {
    TMBFileBED <- data.table::fread(TMBFileBED)
    colnames(TMBFileBED)  <- c("chr", "start", "end")
    assign("TMBFileBED", TMBFileBED, envir = globalenv())
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
  
  # copy number (only 1 per patient with .CNVs extension)
  cnvData <- list.files(path = topDir, pattern = "*.CNVs$", recursive = TRUE, full.names = T)
  if(length(cnvData) == 1){
    cnvData <- data.table::fread(cnvData, header = F)
    cnvData <- as.data.frame(cnvData)
    assign("cnvData", cnvData, envir = globalenv())
  } else {
    # try to read p-value data if available
    cnvData <- list.files(path = topDir, pattern = "*.CNVs.p.value.txt", recursive = TRUE, full.names = T)
    if(length(cnvData) == 1){
      cnvData <- data.table::fread(cnvData, skip = 1, header = F)
      cnvData <- cnvData[,1:5]
      cnvData <- as.data.frame(cnvData)
      assign("cnvData", cnvData, envir = globalenv())
    }
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
  # read fusion files/filter them/merge them
  fusFiles <- list.files(path = topDir, pattern = fusPattern, recursive = TRUE, full.names = T)
  if(length(fusFiles) >= 1){
    fusFiles <- lapply(fusFiles, filterFusions)
    fusData <- do.call('rbind', fusFiles)
    fusData <- unique(fusData)
    if(nrow(fusData) >= 1){
      assign("fusData", fusData, envir = globalenv())
    }
  }
  
  # expression data (can be only 1 per patient with .genes.results)
  expDat <- list.files(path = topDir, pattern = "*.genes.results*", recursive = TRUE, full.names = T)
  if(length(expDat) == 1){
    expData <- read.delim(expDat)
    expData[,"gene_id"] <- sapply(as.character(expData[,1]), FUN = getEns)
    expData <- expData[order(-expData[,"FPKM"]),]
    expData <- expData[!duplicated(expData[,1]),]
    expData[,1] <- sapply(as.character(expData[,1]), FUN = remDotStuff)
    rownames(expData) <- expData[,1]
    assign("expData", expData, envir = globalenv())
  }
}

