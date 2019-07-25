#################################
# High Confidence Alterations P2
#################################

remComma <- function(x) {
  out <- substr(x, nchar(x), nchar(x))
  out <- ifelse(out==",", substr(x, 1, nchar(x)-1), x)
  return(out)
}

getGeneFromMut <- function(x) {
  myGene <-strsplit(x, ":")[[1]][[1]]
  return(myGene)
}

getGeneFromFus <- function(x) {
  myGene1 <-strsplit(x, "_")[[1]][[1]]
  myGene2 <-strsplit(x, "_")[[1]][[2]]
  return(c(myGene1, myGene2))
}

highConfidenceFindingsTable <- function(delRPKM=10) {

  myTable <- allFindingsTable()
  myTable <- myTable[!grepl("Pathway", myTable[,"Type"]),]
  myTable <- myTable[!grepl("Outlier", myTable[,"Type"]),]

  # expression is critical
  if(exists('expData')){
    rnaEvidence <-   RNASeqAnalysisOut[[3]]
    rnaEvidence[,"Gene"] <- rownames(rnaEvidence)
    
    # Get only significant sets
    sigGeneSets <- RNASeqAnalysisOut[[2]][[2]]
    sigGeneSets <- sigGeneSets[sigGeneSets[,"P_VAL"]<0.01,]
    sigGeneSets <- sigGeneSets[,c("Pathway", "Direction")]
    hallMarkSetsTS <- merge(hallMarkSetsTS, sigGeneSets, by.x="ind", by.y="Pathway")
    hallMarkSetsTS[,"ind"] <- paste(hallMarkSetsTS[,"ind"], "(",hallMarkSetsTS[,"Direction"], ")", sep="")
    hallMarkSetsTS <- hallMarkSetsTS[,c("ind", "values")]
    hallMarkSetsTS <- hallMarkSetsTS %>% group_by(values) %>% summarize(ind=paste(ind, collapse=","))
    hallMarkSetsTS <- data.frame(hallMarkSetsTS)
    
    # Supporting Evidence for Deletions
    myTableDel <- myTable[myTable[,2]=="Deletion",]
    if(nrow(myTableDel) > 0){
      myTableDel <- merge(myTableDel, rnaEvidence, by.x="Abberation", by.y="Gene", all.x=T)
      myTableDel <- merge(myTableDel, hallMarkSetsTS, by.x="Abberation", by.y="values", all.x=T)
      myTableDel <- myTableDel[myTableDel[,"SampleX"]<10,]
      if(nrow(myTableDel)>0) {
        myTableDel[,"SupportEv"] <- paste("FPKM=", myTableDel[,"SampleX"], ifelse(is.na(myTableDel[,"ind"]), "", paste(", Pathway: ", myTableDel[,"ind"], sep="")), sep="")
        myTableDel <- myTableDel[,c("Abberation", "Type", "Details", "Score", "Drugs", "Trials", "SupportEv")]
      } else {
        colnames(myTableDel) <- c("Abberation", "Type", "Details", "Score", "Drugs", "Trials", "SupportEv")
      }
    } else {
      myTableDel <- data.frame()
    }
    
    
    # Supporting Evidence for Amplifications
    myTableAmp <- myTable[myTable[,2]=="Amplification",]
    if(nrow(myTableAmp) > 0){
      myTableAmp <- merge(myTableAmp, rnaEvidence, by.x="Abberation", by.y="Gene", all.x=T)
      myTableAmp <- myTableAmp[myTableAmp[,"SampleX"]>100,]
      if(nrow(myTableAmp)>0) {
        myTableAmp[,"SupportEv"] <- paste("FPKM=", myTableAmp[,"SampleX"], ifelse(is.na(myTableAmp[,"ind"]), "", paste(", Pathway: ", myTableAmp[,"ind"], sep="")), sep="")
        myTableAmp <- myTableAmp[,c("Abberation", "Type", "Details", "Score", "Drugs", "Trials", "SupportEv")]
      } else {
        colnames(myTableAmp) <- c("Abberation", "Type", "Details", "Score", "Drugs", "Trials", "SupportEv")
      }
    } else {
      myTableAmp <- data.frame()
    }
    
    
    # Supporting Evidence for Mutations Oncogene - Expression is listed, Pathway is significant
    myTableMut <- myTable[myTable[,2]=="Mutation",]
    if(nrow(myTableMut) > 0){
      myTableMut[,"Gene"] <- sapply(myTableMut[,"Abberation"], FUN=getGeneFromMut)
      myTableMut <- merge(myTableMut, rnaEvidence, by.x="Gene", by.y="Gene", all.x=T)
      myTableMut <- merge(myTableMut, hallMarkSetsTS, by.x="Gene", by.y="values", all.x=T)
      myTableMut[,"SupportEv"] <- paste("FPKM=", myTableMut[,"SampleX"], ifelse(is.na(myTableMut[,"ind"]), "", paste(", Pathway: ", myTableMut[,"ind"], sep="")), sep="")
      myTableMut <- myTableMut[,c("Abberation", "Type", "Details", "Score", "Drugs", "Trials", "SupportEv")]
    } else {
      myTableMut <- data.frame()
    }
    
    
    # Supporting Evidence for Mutations Oncogene - Expression is listed, Pathway is significant
    myTableFus <- myTable[myTable[,2]=="Fusion",]
    if(nrow(myTableFus) > 0){
      myTableFus[,c("Gene1", "Gene2")] <- sapply(myTableFus[,"Abberation"], FUN=getGeneFromFus)
      myTableFus <- merge(myTableFus, rnaEvidence, by.x="Gene1", by.y="Gene", all.x=T)
      myTableFus <- merge(myTableFus, rnaEvidence, by.x="Gene2", by.y="Gene", all.x=T)
      myTableFus <- merge(myTableFus, hallMarkSetsTS, by.x="Gene1", by.y="values", all.x=T)
      myTableFus <- merge(myTableFus, hallMarkSetsTS, by.x="Gene2", by.y="values", all.x=T)
      myTableFus[,"SupportEv"] <- paste("FPKM=", myTableFus[,"SampleX.x"],
                                        ", ",
                                        myTableFus[,"SampleX.y"],
                                        ifelse(is.na(myTableFus[,"ind.x"]), "", paste(", Pathway: ", myTableFus[,"ind.x"], ",", sep="")),
                                        ifelse(is.na(myTableFus[,"ind.y"]), "", paste("", myTableFus[,"ind.y"], sep="")), sep="")
      myTableFus <- myTableFus[,c("Abberation", "Type", "Details", "Score", "Drugs", "Trials", "SupportEv")]
    } else {
      myTableFus <- data.frame()
    }
    
    myTable <- rbind(myTableAmp, myTableDel, myTableMut, myTableFus)
    myTable[,"SupportEv"] <- sapply(myTable[,"SupportEv"], FUN=remComma)
    colnames(myTable)[ncol(myTable)] <- "Supporting Evidence"
    myTable <- unique(myTable)
  } else {
    myTable <- data.frame()
  }
  return(myTable)
}
