#####################
# All Findings Table
#####################

source('code/filterDruggability.R')
source('code/filterFusions.R')
source('code/filterCNV.R')
source('code/filterMutations.R')

allFindingsTable <- function() {
  # Druggability
  drugData <- filterDruggability()
  
  # First Get Mutations
  if(exists('mutData')){
    tmpMut <- filterMutations()
    tmpMut[,"Abberation"] <- ifelse(tmpMut[,"HGVSp_Short"]!="", paste(tmpMut[,"Hugo_Symbol"], tmpMut[,"HGVSp_Short"], sep=": "), as.character(tmpMut[,"Hugo_Symbol"]))
    tmpMut[,"Type"] <- "Mutation"
    tmpMut[,"Details"] <- paste("Mutation Type: ", tmpMut[,"Variant_Classification"], sep="")
    tmpMut[,"Score"] <- "None"
    tmpMut[,"Trials"] <- "None"
    tmpMut <- merge(tmpMut, drugData, by.x="Hugo_Symbol", by.y="gene_name", all.x=T)
    tmpMut <- tmpMut[,c("Abberation", "Type", "Details", "Score", "Drugs", "Trials")]
  } else {
    tmpMut <- data.frame()
  }
  
  # Now Fusions
  if(exists('fusData')){
    if(fusion_method == "arriba"){
      tmpFus <- filterFusions_aribba()
      tmpFus[,"Abberation"] <- gsub("--", "-", tmpFus[,"X.fusion_name"])
      tmpFus[,"Type"] <- "Fusion"
      tmpFus[,"Details"] <- paste("Fusion Type: ",tmpFus[,"Splice_type"], sep="")
      tmpFus[,"Score"] <- "None"
      tmpFus[,"Trials"] <- "None"
      tmpFus <- merge(tmpFus, drugData, by.x="HeadGene", by.y="gene_name", all.x=T)
      tmpFus <- merge(tmpFus, drugData, by.x="TailGene", by.y="gene_name", all.x=T)
      tmpFus[,"Drugs"] <- paste(tmpFus[,"Drugs.x"], tmpFus[,"Drugs.y"], sep=", ")
      tmpFus <- tmpFus[,c("Abberation", "Type", "Details", "Score", "Drugs", "Trials")]
    } else if(fusion_method == "star") {
      tmpFus <- filterFusions_star()
    }
  } else {
    tmpFus <- data.frame()
  }
  
  # Now Copy Number
  if(exists('cnvData')){
    cnvGenes <- createCopyNumber()
    assign("cnvGenes", cnvGenes, envir = globalenv())
    tmpCnv <- filterCNV()
    tmpCnv[,"Abberation"] <- tmpCnv[,1]
    tmpCnv[,"Type"] <- ifelse(tmpCnv[,2]>2, "Amplification", "Deletion")
    tmpCnv[,"Details"] <- paste("Copy Number Value: ",tmpCnv[,2], sep="")
    tmpCnv[,"Score"] <- "None"
    tmpCnv[,"Trials"] <- "None"
    tmpCnv <- merge(tmpCnv, drugData, by.x="Gene", by.y="gene_name", all.x=T)
    tmpCnv <- tmpCnv[,c("Abberation", "Type", "Details", "Score", "Drugs", "Trials")]
  } else {
    tmpCnv <- data.frame()
  }
  
  # Now Expression
  if(exists('RNASeqAnalysisOut')){
    tmpExp <- RNASeqAnalysisOut[[1]][[2]]
    tmpExp[,"Abberation"] <-rownames(tmpExp)
    tmpExp[,"Type"] <- c(rep("Outlier-High (mRNA)", 20), rep("Outlier-Low (mRNA)", 20))
    tmpExp[,"Details"] <- paste("Z-Score / FPKM: ",round(tmpExp[,"Z_Score"],2), " / ", tmpExp[,"FPKM"], sep="")
    tmpExp[,"Score"] <- "None"
    tmpExp[,"Trials"] <- "None"
    tmpExp <- merge(tmpExp, drugData, by.x="Abberation", by.y="gene_name", all.x=T)
    tmpExp <- tmpExp[,c("Abberation", "Type", "Details", "Score", "Drugs", "Trials")]
  } else {
    tmpExp <- data.frame()
  }
  
  # Now Pathway
  if(exists('RNASeqAnalysisOut')){
    tmpPath <- RNASeqAnalysisOut[[2]][[2]]
    tmpPathUp <- tmpPath[tmpPath[,"Direction"]=="Up",][1:20,]
    tmpPathDown <- tmpPath[tmpPath[,"Direction"]=="Down",][1:20,]
    tmpPath <- rbind(tmpPathUp, tmpPathDown)
    tmpPath[,"Abberation"] <-gsub("HALLMARK_", "", tmpPath[,"Pathway"])
    tmpPath[,"Abberation"] <-gsub("_", " ", tmpPath[,"Abberation"])
    tmpPath[,"Type"] <- c(rep("Pathway Up", 20), rep("Pathway Down", 20))
    tmpPath[,"Details"] <- paste("P-Value: ",as.character(formatC(tmpPath[,"P_VAL"], format = "e", digits = 2)), sep="")
    tmpPath[,"Score"] <- "None"
    tmpPath[,"Drugs"] <- "None"
    tmpPath[,"Trials"] <- "None"
    tmpPath <- tmpPath[,c("Abberation", "Type", "Details", "Score", "Drugs", "Trials")]
  } else {
    tmpPath <- data.frame()
  }
  
  # Now Merge Together
  # if(exists('tmpCnv')){
  #   allFindingsDF <- rbind(tmpMut, tmpFus, tmpCnv, tmpExp, tmpPath)
  # } else {
  #   allFindingsDF <- rbind(tmpMut, tmpFus, tmpExp, tmpPath)
  # }
  allFindingsDF <- rbind(tmpMut, tmpFus, tmpCnv, tmpExp, tmpPath)
  if(nrow(allFindingsDF) > 0){
    allFindingsDF[is.na(allFindingsDF[,"Drugs"]),"Drugs"]<- "None"
    allFindingsDF[allFindingsDF[,"Drugs"]=="NA, NA","Drugs"]<- "None"
    allFindingsDF <- unique(allFindingsDF)
  }
  return(allFindingsDF)
}
