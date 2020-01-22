#####################
# All Findings Table
#####################

source('code/filterDruggability.R')
source('code/filterFusions.R')
source('code/filterCNV.R')
source('code/filterMutations.R')
source('code/annotateMutations.R')
source('code/createCopyNumber.R')

allFindingsTable <- function() {
  # Druggability
  drugData <- filterDruggability()
  
  # Somatic Mutations
  if(exists('mutData')){
    tmpMut <- annotateMutations() # don't filter, just annotate
    if(nrow(tmpMut) > 0){
      tmpMut[,"Aberration"] <- ifelse(tmpMut[,"HGVSp_Short"]!="", paste(tmpMut[,"Hugo_Symbol"], tmpMut[,"HGVSp_Short"], sep=": "), as.character(tmpMut[,"Hugo_Symbol"]))
      tmpMut[,"Type"] <- tmpMut[,"Type"]
      tmpMut[,"Details"] <- paste0("Mutation Type: ", tmpMut[,"Variant_Classification"])
      tmpMut[,"Score"] <- "N/A"
      tmpMut[,"Trials"] <- "N/A"
      tmpMut <- merge(tmpMut, drugData, by.x="Hugo_Symbol", by.y="gene_name", all.x=T)
      tmpMut <- tmpMut[,c("Aberration", "Type", "Details", "Score", "Drugs", "Trials")]
    } else {
      tmpMut <- data.frame()
    }
  } else {
    tmpMut <- data.frame()
  }
  
  # Now Fusions
  if(exists('fusData')){
    tmpFus <- fusData
    tmpFus[,"Aberration"] <- tmpFus[,"X.fusion_name"]
    tmpFus[,"Type"] <- "Fusion"
    tmpFus[,"Details"] <- paste0("Fusion Type: ",tmpFus[,"Splice_type"])
    tmpFus[,"Score"] <- "N/A"
    tmpFus[,"Trials"] <- "N/A"
    tmpFus <- merge(tmpFus, drugData, by.x="HeadGene", by.y="gene_name", all.x=T)
    tmpFus <- merge(tmpFus, drugData, by.x="TailGene", by.y="gene_name", all.x=T)
    tmpFus[,"Drugs"] <- paste(tmpFus[,"Drugs.x"], tmpFus[,"Drugs.y"], sep=", ")
    tmpFus <- tmpFus[,c("Aberration", "Type", "Details", "Score", "Drugs", "Trials")]
  } else {
    tmpFus <- data.frame()
  }
  
  # Now Copy Number
  if(exists('cnvData')){
    cnvGenes <- createCopyNumber()
    assign("cnvGenes", cnvGenes, envir = globalenv())
    tmpCnv <- filterCNV()
    if(nrow(tmpCnv) >= 1){
      tmpCnv[,"Aberration"] <- tmpCnv[,1]
      tmpCnv[,"Type"] <- ifelse(tmpCnv[,2]>2, "Amplification", "Deletion")
      tmpCnv[,"Details"] <- paste("Copy Number Value: ",tmpCnv[,2], sep="")
      tmpCnv[,"Score"] <- "N/A"
      tmpCnv[,"Trials"] <- "N/A"
      tmpCnv <- merge(tmpCnv, drugData, by.x="Gene", by.y="gene_name", all.x=T)
      tmpCnv <- tmpCnv[,c("Aberration", "Type", "Details", "Score", "Drugs", "Trials")]
    }  else {
      tmpCnv <- data.frame()
    }
  } else {
    tmpCnv <- data.frame()
  }
  
  # Now Expression
  if(exists('expData')){
    tmpExp <- RNASeqAnalysisOut[[1]][[2]]
    tmpExp[,"Aberration"] <-rownames(tmpExp)
    tmpExp[,"Type"] <- c(rep("Outlier-High (mRNA)", 20), rep("Outlier-Low (mRNA)", 20))
    tmpExp[,"Details"] <- paste("Z-Score / FPKM: ",round(tmpExp[,"Z_Score"],2), " / ", tmpExp[,"FPKM"], sep="")
    tmpExp[,"Score"] <- "N/A"
    tmpExp[,"Trials"] <- "N/A"
    tmpExp <- merge(tmpExp, drugData, by.x="Aberration", by.y="gene_name", all.x=T)
    tmpExp <- tmpExp[,c("Aberration", "Type", "Details", "Score", "Drugs", "Trials")]
  } else {
    tmpExp <- data.frame()
  }
  
  # Now Pathway
  if(exists('expData')){
    tmpPath <- RNASeqAnalysisOut[[2]][[2]]
    tmpPathUp <- tmpPath[tmpPath[,"Direction"]=="Up",][1:20,]
    tmpPathDown <- tmpPath[tmpPath[,"Direction"]=="Down",][1:20,]
    tmpPath <- rbind(tmpPathUp, tmpPathDown)
    tmpPath[,"Aberration"] <-gsub("HALLMARK_", "", tmpPath[,"Pathway"])
    tmpPath[,"Aberration"] <-gsub("_", " ", tmpPath[,"Aberration"])
    tmpPath[,"Type"] <- c(rep("Pathway Up", 20), rep("Pathway Down", 20))
    tmpPath[,"Details"] <- paste("P-Value: ",as.character(formatC(tmpPath[,"P_VAL"], format = "e", digits = 2)), sep="")
    tmpPath[,"Score"] <- "N/A"
    tmpPath[,"Drugs"] <- "N/A"
    tmpPath[,"Trials"] <- "N/A"
    tmpPath <- tmpPath[,c("Aberration", "Type", "Details", "Score", "Drugs", "Trials")]
  } else {
    tmpPath <- data.frame()
  }
  
  # Now Merge Together
  allFindingsDF <- rbind(tmpMut, tmpFus, tmpCnv, tmpExp, tmpPath)
  if(nrow(allFindingsDF) > 0){
    allFindingsDF[is.na(allFindingsDF[,"Drugs"]),"Drugs"]<- "N/A"
    allFindingsDF[allFindingsDF[,"Drugs"]=="NA, NA","Drugs"]<- "N/A"
    allFindingsDF <- unique(allFindingsDF)
  }
  return(allFindingsDF)
}
