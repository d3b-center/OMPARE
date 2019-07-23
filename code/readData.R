#############################
# Read Patient-specific Data
# Save all data to global env
#############################

source('code/createCopyNumber.R')
source('code/helper.R')

readData <- function(topDir, fusion_method = c("star","arriba")){
  
  # patient sample info (at most one with .txt extension)
  sampleInfo <- list.files(path = paste0(topDir, "Clinical"), pattern = "*.txt", full.names = T)
  if(length(sampleInfo) == 1){
    sampleInfo <- read.delim(sampleInfo, stringsAsFactors = F)
  } else {
    # assign n/a if no clinical info is present
    sampleInfo <- setNames(data.frame(t(rep("n/a", 15)), stringsAsFactors = F), 
                           nm = c("patientName","reportDate","reportVersion","primRelapse",
                                  "tumorType","tumorLocation","collectionDate","labDirector",
                                  "pathologist","primPhysician","medicalFacility","ethnicity",
                                  "dob","sex","subjectID"))
  }
  assign("sampleInfo", sampleInfo, envir = globalenv())
  
  # mutation data (can be multiple with .maf extension)
  mutFiles <- list.files(path = topDir, pattern = "*.maf", recursive = TRUE, full.names = T)
  if(length(mutFiles) >= 1){
    mutFiles <- lapply(mutFiles, data.table::fread, skip = 1, stringsAsFactors = F)
    mutData <- data.table::rbindlist(mutFiles)
    mutData <- as.data.frame(mutData)
    mutData <- unique(mutData)
    assign("mutData", mutData, envir = globalenv())
  } 
  
  # copy number (only 1 per patient with .CNVs extension)
  cnvData <- list.files(path = topDir, pattern = "*.CNVs", recursive = TRUE, full.names = T)
  if(length(cnvData) == 1){
    cnvData <- data.table::fread(cnvData, header = F)
    cnvData <- as.data.frame(cnvData)
    assign("cnvData", cnvData, envir = globalenv())
  }
  
  # fusion data (currently either star or arriba implemented)
  if(fusion_method == "star"){
    fusPattern = "*star-fusion.fusion_candidates.final"
  } else if(fusion_method == "arriba") {
    fusPattern = "*.arriba.fusions.tsv"
  }
  fusFiles <- list.files(path = topDir, pattern = fusPattern, recursive = TRUE, full.names = T)
  if(length(fusFiles) >= 1){
    fusFiles <- lapply(fusFiles, read.delim)
    fusData <- do.call('rbind', fusFiles)
    fusData <- unique(fusData)
    assign("fusData", fusData, envir = globalenv())
  }
  
  # expression data (can be only 1 per patient with .genes.results)
  expDat <- list.files(path = topDir, pattern = "*.genes.results", recursive = TRUE, full.names = T)
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

