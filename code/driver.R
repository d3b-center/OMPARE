###########################################
#Purpose: Main code to generate all tables/figures for report
#Author: Pichai Raman
#Date: 3/21/2019
###########################################

#############################
#Load packages & source-
#############################
library("tidyverse", verbose=F, quietly=T)
library("copynumber", verbose=F, quietly=T)
library("RCircos", verbose=F, quietly=T)
library("biomaRt", verbose=F, quietly=T)
library("GeneNetworkBuilder", verbose=F, quietly=T)
library("DT", verbose=F, quietly=T)
library("GSVA", verbose=F, quietly=T);
library("GSEABase", verbose=F, quietly=T);
library("igraph", verbose=F, quietly=T);
library("network", verbose=F, quietly=T);
library("sna", verbose=F, quietly=T);
library("ggnetwork", verbose=F, quietly=T);
library("Rtsne", verbose=F, quietly=T);
library("RColorBrewer", verbose=F, quietly=T);
library("randomcoloR", verbose=F, quietly=T);
library("survminer", verbose=F, quietly=T);
library("RTCGA.clinical", verbose=F, quietly=T);
library("survival", verbose=F, quietly=T);
library("optparse", verbose=F, quietly=T);
library("xCell", verbose=F, quietly=T);
library("decompTumor2Sig", verbose=F, quietly=T);
library("stringr", verbose=F, quietly=T);

source("code/dgidbParse.R")
source("code/parseSurvival.R");
source("code/copyNumberConvert.R")
source("code/patientSampleInfo.R")
source("code/pubTheme.R");
source("code/germlineAnalysis.R")
#############################
#-End Load packages & source
#############################


#############################
#-Read Data
#############################
#Mutation Data
mutData <- read.delim("../data/PNOC008/PNOC008-2/f39e47a9-8a9d-4ebb-b8ec-40569961f7b2.strelka2_somatic.vep.maf", skip=1, stringsAsFactors = F);
mutData2 <- read.delim("../data/PNOC008/PNOC008-2/f39e47a9-8a9d-4ebb-b8ec-40569961f7b2.mutect2_somatic.vep.maf", skip=1, stringsAsFactors = F);
mutData <- rbind(mutData, mutData2[,colnames(mutData)])
mutData <- unique(mutData);

#Copy Number Data
#cnvData <- read.delim("../data/CBTTC-HGG/CNV/ae4ce725-d3c4-455d-822e-6c5067444b5e.CNVs", header=F);
#cnvGenes <- createCopyNumber()




#Expression Data
expData <- read.delim("../data/PNOC008/PNOC008-2/eb06d52a-110b-4c0e-9e90-94ea68d2f698.rsem.genes.results")
getEns <- function(x)
{
  out <- strsplit(x, split="_")[[1]][1]
  return(out);
}
expData[,"gene_id"] <- sapply(as.character(expData[,1]), FUN=getEns)
expData <- expData[order(-expData[,"FPKM"]),]
remDotStuff <- function(x)
{
  out <- strsplit(x, "\\.")[[1]][1]
}
expData <- expData[!duplicated(expData[,1]),]
expData[,1] <- sapply(as.character(expData[,1]), FUN=remDotStuff)
rownames(expData) <- expData[,1];

source("code/RNASeqAnalysis.R")


#Fusion Data
#fusData <- read.delim("../data/PNOC008/PNOC008-2/7316-325.local.transcript.converted.pe.star-fusion.fusion_candidates.final")

fusData <- read.delim("../data/PNOC008/PNOC008-2/eb06d52a-110b-4c0e-9e90-94ea68d2f698.arriba.fusions.tsv")

##Reference Data
cancerGenes <- read.delim("../data/Reference/CancerGeneList.tsv")
tsgGenes <- read.delim("../data/Reference/Human_TSGs.txt")
chrMap <- read.delim("../data/Reference/mart_export_genechr_mapping.txt", stringsAsFactors =F)
geneMania <- read.delim("../data/Reference/GeneManiaNetwork.txt", stringsAsFactors =F)
diseaseSpecificFields <- read.delim("../data/Reference/DiseaseSpecificFields.txt")
hallMarkSets <- getGmt("../data/Reference/mSigDB/h.all.v6.2.symbols.gmt", collectionType=BroadCollection(), geneIdType= SymbolIdentifier())
hallMarkSets <- geneIds(hallMarkSets);
hallMarkSetsTS <- stack(hallMarkSets)
load("../data/Reference/cbttc_genes_fpkm_1110.RData");
res[,1] <- sapply(as.character(res[,1]), FUN=remDotStuff)
mapping <- read.delim("../data/Reference/mappingFile.txt", header=F, stringsAsFactors=F);
clinData <- read.delim("../data/Reference/study_view_clinical_data.txt", stringsAsFactors=F);
tmbData <- read.csv("../data/reference/complete_results.csv", stringsAsFactors=F);
parseSurvival();
survData <- read.delim("../data/Reference/survData.txt", stringsAsFactors=F);

#############################
#-End Read Data
#############################


####################################################
#Function to filter mutations-
####################################################
filterMutations <- function(myMutData=mutData, myCancerGenes=cancerGenes)
{
  mutDataFilt <- myMutData;
  
  #Filter to only protein coding genes
  mutDataFilt <- mutDataFilt[mutDataFilt[,"BIOTYPE"]=="protein_coding",]
  
  #Filter by Variant Classification
  keepVC <- c("Missense_Mutation", "Splice_Region", "3'UTR", "5'UTR", "In_Frame_Del")
  mutDataFilt <- mutDataFilt[mutDataFilt[,"Variant_Classification"]%in%keepVC,]

  #Filter by Variant IMPACT
  keepVI <- c("MODIFIER", "MODERATE", "HIGH")
  mutDataFilt <- mutDataFilt[mutDataFilt[,"IMPACT"]%in%c(keepVI),]

  #Filter by Cancer Gene Census
  myCancerGenes <- as.character(myCancerGenes[,1]);
  mutDataFilt <- mutDataFilt[mutDataFilt[,"Hugo_Symbol"]%in%myCancerGenes,];
  return(mutDataFilt);
}
####################################################
#-End Function to filter mutations
####################################################


####################################################
#Function to filter fusions-
####################################################
filterFusions_star <- function(myFusData=fusData, myCancerGenes=cancerGenes, myJunctionReads=2)
{
  fusDataFilt <- myFusData;
  
  #FunctionToSplitGenes
  splitMyGene <- function(x)
  {
    x <- strsplit(x, split="\\^")[[1]][1]
    return(x);
  }

  #Filter by Number of Reads
  fusDataFilt <- fusDataFilt[fusDataFilt[,"JunctionReads"]>myJunctionReads,]
  
  #Filter by Cancer Gene Census
  myCancerGenes <- as.character(myCancerGenes[,1]);
  fusDataFilt[,"HeadGene"] <- sapply(as.character(fusDataFilt[,"LeftGene"]), FUN=splitMyGene)
  fusDataFilt[,"TailGene"] <- sapply(as.character(fusDataFilt[,"RightGene"]), FUN=splitMyGene)
  fusDataFilt <- fusDataFilt[fusDataFilt[,"HeadGene"]%in%myCancerGenes|fusDataFilt[,"TailGene"]%in%myCancerGenes,];
  
  return(fusDataFilt);
  
}

filterFusions_aribba <- function(myFusData=fusData, myCancerGenes=cancerGenes)
{
  fusDataFilt <- myFusData;
  fusDataFilt[,"X.fusion_name"] <- paste(as.character(fusDataFilt[,"X.gene1"]), "-", as.character(fusDataFilt[,"gene2"]), sep="");
  fusDataFilt[,"Splice_type"] <- fusDataFilt[,"type"]
  
  
  #Filter by Number of Reads
  fusDataFilt <- fusDataFilt[fusDataFilt[,"confidence"]!="low",]
  
  #Filter by Cancer Gene Census
  myCancerGenes <- as.character(myCancerGenes[,1]);
  fusDataFilt[,"HeadGene"] <- as.character(fusDataFilt[,"X.gene1"])
  fusDataFilt[,"TailGene"] <- as.character(fusDataFilt[,"gene2"])
  fusDataFilt <- fusDataFilt[fusDataFilt[,"HeadGene"]%in%myCancerGenes|fusDataFilt[,"TailGene"]%in%myCancerGenes,];
  
  return(fusDataFilt);
  
}

####################################################
#-End Function to filter fusions
####################################################

####################################################
#Function to filter CNV-
####################################################
filterCNV <- function(myCNVData=cnvGenes, myCancerGenes=cancerGenes, myTSGenes=tsgGenes, cutoffHigh=4, cutoffLow=1)
{
  cnvDataFilt <- myCNVData;
 
  #Filter by Cancer Gene Census
  myTSGenes <- as.character(myTSGenes[,2]);
  myOncogenes <- setdiff(as.character(myCancerGenes[,1]), myTSGenes);
  
  cnvDataFiltUp <- cnvDataFilt[cnvDataFilt[,2]>cutoffHigh,]
  cnvDataFiltUp <- cnvDataFiltUp[cnvDataFiltUp[,1]%in%myOncogenes,]

  cnvDataFiltDown <- cnvDataFilt[cnvDataFilt[,2]<cutoffLow,]
  cnvDataFiltDown <- cnvDataFiltDown[cnvDataFiltDown[,1]%in%myTSGenes,]
  cnvDataFilt <- rbind(cnvDataFiltUp, cnvDataFiltDown);

  return(cnvDataFilt);
  
}
####################################################
#-End Function to filter CNV
####################################################



####################################################
####################################################
#Key Clinical Findings P1
####################################################
####################################################

  
patientSampleInfo <- function()
{
  df1 <- data.frame(c("Subject ID", "Sex", "DOB", "Ethnicity"), c(getSubjectID(), getSex(), getDOB(), getEthnicity()))
  df2 <- data.frame(c("Medical Facility", "Primary Physician", "Pathologist", "Lab Director"), c(getMedicalFacility(), getPrimPhysician(), getPathologist(), getLabDirector()))
  df3 <- data.frame(c("Collection Date", "Tumor Location", "Tumor Type", "P/R"), c(getCollectionDate(), getTumorLocation(), getTumorType(), getPrimRelapse()))
  return(cbind(df1,df2, df3));
}

keyClinicalFindingsTable <- function()
{
  return(highConfidenceFindingsTable())
}

diseaseSpecificInformation <- function()
{
  tmpGeneFindings <- allFindingsTable();
  tmpGeneFindings <- tmpGeneFindings[!grepl("Pathway", tmpGeneFindings[,"Type"]),]
  diseaseSpecificFields
  
  #Check everything 
  getStatus <- function(x)
  {
    tmpGenes <- x[[3]]
    tmpGenes <- trimws(strsplit(tmpGenes, ",")[[1]])
    tmpOut <- sapply(tmpGenes, FUN=grepl, x=tmpGeneFindings[,1])
    paste(paste(tmpGeneFindings[as.logical(rowSums(tmpOut)),1], ":", tmpGeneFindings[as.logical(rowSums(tmpOut)),2], sep=""), collapse=", ");
  }
  diseaseSpecificFields[,"Value"] <- apply(diseaseSpecificFields, FUN=getStatus, MARGIN=1)
  diseaseSpecificFields[diseaseSpecificFields[,"Value"]==":","Value"] <- "Normal";
  diseaseSpecificFields <- diseaseSpecificFields[,c("Field_name", "Value")]
  return(diseaseSpecificFields)
}


germlineInformation <- function()
{
  tmpGeneFindings <- allFindingsTable();
  tmpGeneFindings <- tmpGeneFindings[!grepl("Pathway", tmpGeneFindings[,"Type"]),]
  diseaseSpecificFields
  
  #Check everything 
  getStatus <- function(x)
  {
    tmpGenes <- x[[3]]
    tmpGenes <- trimws(strsplit(tmpGenes, ",")[[1]])
    tmpOut <- sapply(tmpGenes, FUN=grepl, x=tmpGeneFindings[,1])
    paste(paste(tmpGeneFindings[as.logical(rowSums(tmpOut)),1], ":", tmpGeneFindings[as.logical(rowSums(tmpOut)),2], sep=""), collapse=", ");
  }
  diseaseSpecificFields[,"Value"] <- apply(diseaseSpecificFields, FUN=getStatus, MARGIN=1)
  diseaseSpecificFields[diseaseSpecificFields[,"Value"]==":","Value"] <- "Normal";
  diseaseSpecificFields <- diseaseSpecificFields[,c("Field_name", "Value")]
  return(diseaseSpecificFields)
}



genomicSummary <- function()
{
  tmpRightHead <- c("High Confidence Genomic Alterations", "Total Genomic Alterations", "Transcriptomic Alterations", "Proteomic Alterations", "Aberrant Pathway Activity")
  
  numLesions <- allFindingsTable();
  numLesions <- nrow(numLesions[numLesions[,"Type"]%in%c("Mutation", "Fusion", "Amplification", "Deletion"),]);
  numTranscripts <- RNASeqAnalysisOut[[1]][[2]];
  numTranscripts <- nrow(numTranscripts[numTranscripts[,1]>3,])
  numProeins <- "NA"
  numPathways <- RNASeqAnalysisOut[[2]][[2]];
  numPathways <- numPathways[numPathways[,"P_VAL"]<0.01,]
  numPathways <- nrow(numPathways);
 
  tmpVals <- c(nrow(highConfidenceFindingsTable()),numLesions, numTranscripts, numProeins, numPathways);
  df1 <- data.frame(tmpRightHead, tmpVals);
  return(df1);
}


####################################################
####################################################
#End Key Clinical Findings P1
####################################################
####################################################




####################################################
####################################################
#High Confidence Alterations P2
####################################################
####################################################


highConfidenceFindingsTable <- function(delRPKM=10)
{
  
  myTable <- allFindingsTable();
  myTable <- myTable[!grepl("Pathway", myTable[,"Type"]),]
  myTable <- myTable[!grepl("Outlier", myTable[,"Type"]),]

  rnaEvidence <-   RNASeqAnalysisOut[[3]];
  rnaEvidence[,"Gene"] <- rownames(rnaEvidence);
  
  #Get only significant sets
  sigGeneSets <- RNASeqAnalysisOut[[2]][[2]];
  sigGeneSets <- sigGeneSets[sigGeneSets[,"P_VAL"]<0.01,]
  sigGeneSets <- sigGeneSets[,c("Pathway", "Direction")];
  hallMarkSetsTS <- merge(hallMarkSetsTS, sigGeneSets, by.x="ind", by.y="Pathway")
  hallMarkSetsTS[,"ind"] <- paste(hallMarkSetsTS[,"ind"], "(",hallMarkSetsTS[,"Direction"], ")", sep="");
  hallMarkSetsTS <- hallMarkSetsTS[,c("ind", "values")];
  hallMarkSetsTS <- hallMarkSetsTS %>% group_by(values) %>% summarize(ind=paste(ind, collapse=","))
  hallMarkSetsTS <- data.frame(hallMarkSetsTS);
  #Supporting Evidence for Deletions
  myTableDel <- myTable[myTable[,2]=="Deletion",]
  myTableDel <- merge(myTableDel, rnaEvidence, by.x="Abberation", by.y="Gene", all.x=T)
  myTableDel <- merge(myTableDel, hallMarkSetsTS, by.x="Abberation", by.y="values", all.x=T)
  myTableDel <- myTableDel[myTableDel[,"SampleX"]<10,]
  if(nrow(myTableDel)>0)
  {
    myTableDel[,"SupportEv"] <- paste("FPKM=", myTableDel[,"SampleX"], ifelse(is.na(myTableDel[,"ind"]), "", paste(", Pathway: ", myTableDel[,"ind"], sep="")), sep="")
    myTableDel <- myTableDel[,c("Abberation", "Type", "Details", "Score", "Drugs", "Trials", "SupportEv")]
  }
  if(nrow(myTableDel)==0)
  {
    colnames(myTableDel) <- c("Abberation", "Type", "Details", "Score", "Drugs", "Trials", "SupportEv")
  }
  
  #Supporting Evidence for Amplifications
  myTableAmp <- myTable[myTable[,2]=="Amplification",]
  myTableAmp <- merge(myTableAmp, rnaEvidence, by.x="Abberation", by.y="Gene", all.x=T)
  myTableAmp <- myTableAmp[myTableAmp[,"SampleX"]>100,]
  if(nrow(myTableAmp)>0)
  {
    myTableAmp[,"SupportEv"] <- paste("FPKM=", myTableAmp[,"SampleX"], ifelse(is.na(myTableAmp[,"ind"]), "", paste(", Pathway: ", myTableAmp[,"ind"], sep="")), sep="")
    myTableAmp <- myTableAmp[,c("Abberation", "Type", "Details", "Score", "Drugs", "Trials", "SupportEv")]
  }
  if(nrow(myTableAmp)==0)
  {
    colnames(myTableAmp) <- c("Abberation", "Type", "Details", "Score", "Drugs", "Trials", "SupportEv")
  }
  
  #Supporting Evidence for Mutations Oncogene - Expression is listed, Pathway is significant
  myTableMut <- myTable[myTable[,2]=="Mutation",]
  getGeneFromMut <- function(x)
  {
    myGene <-strsplit(x, ":")[[1]][[1]]
    return(myGene);
  }
  myTableMut[,"Gene"] <- sapply(myTableMut[,"Abberation"], FUN=getGeneFromMut);
  myTableMut <- merge(myTableMut, rnaEvidence, by.x="Gene", by.y="Gene", all.x=T)
  myTableMut <- merge(myTableMut, hallMarkSetsTS, by.x="Gene", by.y="values", all.x=T)
  myTableMut[,"SupportEv"] <- paste("FPKM=", myTableMut[,"SampleX"], ifelse(is.na(myTableMut[,"ind"]), "", paste(", Pathway: ", myTableMut[,"ind"], sep="")), sep="")
  myTableMut <- myTableMut[,c("Abberation", "Type", "Details", "Score", "Drugs", "Trials", "SupportEv")]

  #Supporting Evidence for Mutations Oncogene - Expression is listed, Pathway is significant
  myTableFus <- myTable[myTable[,2]=="Fusion",]
  getGeneFromFus <- function(x)
  {
    myGene1 <-strsplit(x, "-")[[1]][[1]]
    myGene2 <-strsplit(x, "-")[[1]][[2]]
    return(c(myGene1, myGene2));
  }
  myTableFus[,c("Gene1", "Gene2")] <- sapply(myTableFus[,"Abberation"], FUN=getGeneFromFus);
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

  remComma <- function(x)
  {
    out <- substr(x, nchar(x), nchar(x));
    out <- ifelse(out==",", substr(x, 1, nchar(x)-1), x);
    return(out);
  }
# myTable <- rbind(myTableAmp, myTableDel, myTableMut, myTableFus)
  myTable <- rbind(myTableMut, myTableFus)
  myTable[,"SupportEv"] <- sapply(myTable[,"SupportEv"], FUN=remComma);
  
  colnames(myTable)[ncol(myTable)] <- "Supporting Evidence";
  myTable <- unique(myTable);
  return(myTable);  
}

####################################################
#Function to get RNA-Seq and Pathway Analysis
####################################################
RNASeqAnalysisOut <- runRNASeqAnalysis(expData)
plotGenes <- function(myRNASeqAnalysisOut=RNASeqAnalysisOut)
{
  geneData <- myRNASeqAnalysisOut[[1]][[2]]
  geneData[,"Direction"] <- ifelse(geneData[,"Z_Score"]>0, "Up", "Down");
  geneData[,"Gene"] <- rownames(geneData);
  geneData <- geneData[order(geneData[,"Z_Score"]),]
  geneData[,"Gene"] <- factor(geneData[,"Gene"], levels=geneData[,"Gene"])
  p <- ggplot(geneData, aes(factor(Gene), y=Z_Score, fill=Direction))+geom_bar(stat="identity")+coord_flip()+theme_bw();
  p <- p+xlab("Gene Symbol")+scale_fill_manual(values = c("forest green", "red"));
  return(p);
}
plotPathway <- function(myRNASeqAnalysisOut=RNASeqAnalysisOut)
{
  pathData <- myRNASeqAnalysisOut[[2]][[2]];
  pathDataUp <- pathData[pathData[,"Direction"]=="Up",]
  pathDataUp[,"Pathway"] <- rownames(pathDataUp);
  pathDataUp <- pathDataUp[order(pathDataUp[,"P_VAL"]),]
  pathDataUp[,"Pathway"] <- factor(pathDataUp[,"Pathway"], levels=pathDataUp[,"Pathway"])
  pathDataUp <- pathDataUp[1:10,];

  pathDataDown <- pathData[pathData[,"Direction"]=="Down",]
  pathDataDown[,"Pathway"] <- rownames(pathDataDown);
  pathDataDown <- pathDataDown[order(pathDataDown[,"P_VAL"]),]
  pathDataDown[,"Pathway"] <- factor(pathDataDown[,"Pathway"], levels=pathDataDown[,"Pathway"])
  pathDataDown <- pathDataDown[1:10,];
  pathData <- rbind(pathDataDown, pathDataUp)
  pathData[,"Direction"] <- factor(pathData[,"Direction"], levels=c("Down", "Up"))
  p <- ggplot(pathData, aes(factor(Pathway), y=(-1)*log10(P_VAL), fill=Direction))+geom_bar(stat="identity")+coord_flip()+theme_bw();
  p <- p+xlab("Pathway Name")+scale_fill_manual(values = c("forest green", "red"))+ylab("-log10 P-Value");
  return(p);
}
####################################################
#End Function to get RNA-Seq and Pathway Analysis
####################################################

####################################################
####################################################
#End High Confidence Alterations P2
####################################################
####################################################







####################################################
####################################################
#Immune Signatures P3
####################################################
####################################################

#Immune Profile
ImmuneProfile <- function()
{
#  resAll <- res;
#  rownames(resAll) <- resAll[,1];
#  combGenes <- intersect(rownames(resAll), rownames(expData))
#  resAll <- cbind(resAll[combGenes,], expData[combGenes,"FPKM"]);
#  colnames(resAll)[ncol(resAll)] <- "PatSample";
#  resAll[,"max"] <- apply(resAll[3:ncol(resAll)], FUN=max, MARGIN=1)
#  resAll <- resAll[order(-resAll[,"max"]),]
#  resAll <- resAll[!duplicated(resAll[,2]),]
#  rownames(resAll) <- resAll[,2];
#  resAll <- resAll[-1:-2];
#  resAll <- resAll[-ncol(resAll)]; 
  

#raw.scores = rawEnrichmentAnalysis(as.matrix(resAll),
#                                   xCell.data$signatures,
#                                   xCell.data$genes)
#  raw.scores[,"CellType"] <- rownames(raw.scores)
#  write.table(raw.scores, "rawScores.txt", sep="\t", row.names=F)

  raw.scores <- read.delim("rawScores.txt")
  raw.scoresTS <- gather(raw.scores, "Sample", "Score", -CellType)
  raw.scoresTS[,"IsSample"] <- ifelse(grepl("PatSample", raw.scoresTS[,"Sample"]), T, F);
  p <- ggplot(raw.scoresTS, aes(CellType, Score))+geom_boxplot(outlier.shape=NA)+theme_bw()+ theme(axis.text.x = element_text(angle = 75, hjust = 1));
  raw.scoresTSSample <- raw.scoresTS[raw.scoresTS[,4]==T,]
  p <-p+geom_point(data=raw.scoresTSSample, aes(CellType, Score), colour = "red", size=3, shape="triangle");
  p <- p+theme(axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14,face="bold"))
  return(p);
}


####################################################
####################################################
#End Immune/Tumor Signatures P3
####################################################
####################################################


####################################################
####################################################
#TMB Tumor Signatures P4
####################################################
####################################################

#TMB Profile
tmbProfile <- function(myTMB=tmbData)
{
  
  #Determine if sample is an outlier in dataset
  isOutlier <- function(distTmp)
  {
    c25 <- as.numeric(quantile(distTmp, .25));
    c75 <- as.numeric(quantile(distTmp, .75));
    myThresh <- (c75-c25)+c75;
    out <- distTmp>myThresh;
  }
  
  
  getMMRFreq <- function(tmpAnnot)
  {
    #Define the Hypermutation profile
    tmpAnnot[,"Outlier"] <- isOutlier(tmpAnnot[,"MUTATION_COUNT"])
    tmpAnnot[,"EST_MUT_PER_MB"] <- tmpAnnot[,"MUTATION_COUNT"]/30;
    return(tmpAnnot);
  }
  number_ticks <- function(n) {function(limits) pretty(limits, n)}
  
  
  ######################
  ##Read data 
  ####################
  
  #Read in data 26504 to start
  
  ######################
  ##Filter data 
  ####################
  
  #Filter to studies that have mutation counts 15373
  myTMB[,"MUTATION_COUNT"] <- as.numeric(myTMB[,"MUTATION_COUNT"]);
  myTMB <- myTMB[!is.na(myTMB[,"MUTATION_COUNT"]),]
  myTMB[,"NAME"] <- as.character(myTMB[,"NAME"])
  
  #Remove PPTP and FMI
  myTMB <- filter(myTMB, DISEASE.NAME!="Pediatric Preclinical Testing Program" & DISEASE.NAME!="Mixed Cancer Types")
  
  #Add Whether Pediatric or Adult
  myTMB[,"Type"] <- "Adult"
  pediatricCancers <- c("Neuroblastoma", "Medulloblastoma", "Acute Myeloid Leukemia", "Pediatric Ganglio", "Diffuse Intrinsic Pontine Glioma")
  myTMB[myTMB[,"DISEASE.NAME"]%in%pediatricCancers, "Type"] <- "Pediatric";
  
  #Sort and order it in decreasing fashion 
  grouped <- group_by(myTMB, Type, DISEASE.NAME)
  tmpOut <- summarise(grouped, median=median(MUTATION_COUNT))
  tmpOut <- arrange(tmpOut, Type, desc(median))
  myOrder <- as.character(as.data.frame(tmpOut)[,"DISEASE.NAME"])
  myTMB[,"DISEASE.NAME"] <- factor(myTMB[,"DISEASE.NAME"], levels=myOrder);
  myTMB[,"EST_MUT_PER_MB"] <- myTMB[,"MUTATION_COUNT"]/30;
  
  #Filter mutations to get TMB of samples
  filtMut <- mutData 
  filtMut <- filtMut[filtMut[,"HGVSp_Short"]!="",];
  tmbSample <- nrow(filtMut)/30;
  
  #Plot it
  p <- ggplot(myTMB, aes(DISEASE.NAME, EST_MUT_PER_MB, fill=Type))+geom_boxplot()+theme_bw();
  p <- p+scale_y_log10(breaks=c(.25, 1, 10, 100, 500))+scale_fill_manual(values=c("blue", "red"))
  p <- p+xlab("Disease")+ylab("Est Mutation Count per MB (Log Scale)")+theme(axis.text.x = element_text(angle = -90, hjust = (0)));
  p <- p+ggtitle("Mutation Count versus Disease")+geom_hline(yintercept=tmbSample, linetype=2)
  return(p)
  
}

tumorSignaturePlot <- function(x)
{
  
  #Have to write out to read in?
  mpfData <- mutData[, c("Chromosome", "Start_Position", "Match_Norm_Seq_Allele1", "Tumor_Seq_Allele2")]
  mpfData <- data.frame("Patient", mpfData);
  write.table(mpfData, "../data/mpfDataFormat.txt", sep="\t", row.names=F, quote=F, col.names=F)
  
  # load the reference genome and the transcript annotation database
  refGenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  transcriptAnno <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  
  # Read in the genome
  genomes <- readGenomesFromMPF("../data/mpfDataFormat.txt", numBases=3, type="Alexandrov",
                                trDir=F, refGenome=refGenome, verbose=F)
  
  signatures <- readAlexandrovSignatures("../data/Reference/signatures_probabilities.txt");
  exposure <- decomposeTumorGenomes(genomes, signatures)[[1]]
  
  exposure <- data.frame(names(exposure), exposure);
  colnames(exposure) <- c("Signature", "Value");
  exposure[,1] <- gsub("\\.", "-", exposure[,1])
  exposure[,1] <- factor(exposure[,1], levels=exposure[,1]);
  p <- ggplot(exposure, aes(Signature, Value))+geom_bar(stat="identity")+theme_bw()+coord_flip();
  p <- p+xlab("Mutational Signature")+ylab("Exposures (percent contribution)")
  return(p)
  
}


####################################################
####################################################
#End TMB Tumor Signatures P4
####################################################
####################################################



####################################################
####################################################
#Genomically Similar Samples P5
####################################################
####################################################

clinData <- merge(clinData, mapping, by.x="Sample.ID", by.y="V1")
clinDataOrig <- clinData;
clinData <- clinData[!grepl("CL", clinData[,"CBTTC_PAIRED_IDS"]),]
clinData <- clinData[,c("V2", "Cancer.Type", "Cancer.Type.Detailed", "TUMOR_TISSUE_SITE")]
clinData <- unique(clinData);
clinData <- clinData[!duplicated(clinData[,"V2"]),]
rownames(clinData)<- clinData[,"V2"]
patVector <- c("PatSample", "High-grade glioma/astrocytoma (WHO grade III/IV)", "High-grade glioma/astrocytoma (WHO grade III/IV)", "Unknown")
clinData <- rbind(clinData, patVector)
rownames(clinData)[nrow(clinData)] <- "PatSample";

##########################
#Filter Expression Data
##########################
#Filter data and get it by gene, remove all genes with "-", ".", "_"

rownames(res) <- res[,1];
combGenes <- intersect(rownames(res), rownames(expData))
res <- cbind(res[combGenes,], expData[combGenes,"FPKM"]);
colnames(res)[ncol(res)] <- "PatSample";
res[,2] <- as.character(res[,2]);
res <- res[!grepl("-", res[,2]),]
res <- res[!grepl("\\.", res[,2]),]
res <- res[!grepl("_", res[,2]),]
resAll <- res;

#Now remove genes that have less 20 FPKM
maxVals <- apply(res[3:ncol(res)], FUN=max, MARGIN=1)
res <- res[maxVals>50,]

#Now take gene with max value
res[,"max"] <- apply(res[3:ncol(res)], FUN=max, MARGIN=1)
res <- res[order(-res[,"max"]),]
res <- res[!duplicated(res[,2]),]
rownames(res) <- res[,2];
res <- res[-1:-2];
res <- res[-ncol(res)]; #Remove max
res <- res[,rownames(clinData)]

##########################
#Get Annotation data ready and constrain columns
##########################

#Get most variable genes
myCV <- function(x) { sd(x)/mean(x)}
myCVs <- apply(res, FUN=myCV, MARGIN=1)
resTmp <- res
resTmp["CV"] <- myCVs;
resTmp <- resTmp[order(-resTmp[,"CV"]),];
resTmp <- resTmp[1:10000,]
resTmp <- resTmp[-ncol(resTmp)]; #Remove cv

getTSNEPlot <- function()
{
 
  tsneOut <- Rtsne(t(log2(resTmp+1)), initial_dims=100, perplexity=30, check_duplicates = FALSE)
  tsneData <- data.frame(tsneOut$Y, colnames(resTmp))
  tsneData <- cbind(clinData, tsneData);
  tsneData[,"SampleX"] <- ifelse(tsneData[,"V2"]=="PatSample", 2, 1);
  tmpCol <- as.character(clinData[,"Cancer.Type"])
  p <- ggplot(tsneData, aes(X1, X2, color=Cancer.Type, size=SampleX, shape=as.character(SampleX)))+geom_jitter(width = 0.5, height = 0.5)+theme_bw()+ggtitle("T-SNE PBTA RNA-Sequencing");
  p <- p+theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
               legend.title=element_text(size=12), 
               legend.text=element_text(size=12))+ guides(shape=FALSE, size=FALSE)
  p;  
}

clinDataHGG <- clinData[grepl("High-grade glioma", clinData[,2]),];
resTmpHGG <- resTmp[,rownames(clinDataHGG)]
allCor <- cor(resTmpHGG[ncol(resTmpHGG)], resTmpHGG)
allCor <- data.frame(t(allCor))
allCor[,"samps"] <- rownames(allCor);
allCor <- allCor[!grepl("Pat", rownames(allCor)),]
allCor <- allCor[intersect(rownames(allCor), survData[,"samps"]),]
allCor <- allCor[order(-allCor[,1]),]

getKMPlot <- function(numNeighbors=15)
{
  mySamps <- allCor[1:numNeighbors,"samps"]
  survData[,"group"] <- survData[,"samps"]%in%mySamps
  survData[,"group"] <- ifelse(survData[,"group"], "Cluster With Patient", "Cluster away from Patient")
  survData <- survData[survData[,1]%in%rownames(allCor),]
  fit <- survfit(Surv(OS_Time, OS_Event) ~ group, data = survData)
  ggsurvplot(
    fit,                     # survfit object with calculated statistics.
    data = survData,  # data used to fit survival curves. 
    risk.table = TRUE,       # show risk table.
    pval = TRUE,             # show p-value of log-rank test.
    conf.int = F,         # show confidence intervals for 
    # point estimaes of survival curves.
    xlim = c(0,1000),        # present narrower X axis, but not affect
    # survival estimates.
    break.time.by = 100,     # break X axis in time intervals by 500.
    ggtheme = theme_minimal(), # customize plot and risk table with a theme.
    risk.table.y.text.col = T, # colour risk table text annotations.
    risk.table.y.text = FALSE # show bars instead of names in text annotations
    # in legend of risk table
  )  
  
}

getSimilarPatients <- function(numNeighbors=15)
{
  mySamps <- allCor[1:numNeighbors,"samps"]
  clinDataTmp <- clinDataOrig[clinDataOrig[,"V2"]%in%mySamps,]
  clinDataTmp <- clinDataTmp[,c("Sample.ID", "Cancer.Type", "AGE", "ETHNICITY", "RACE", "TUMOR_SITE")]
  clinDataTmp[,"Report"] <- paste("<a href='www.google.com'>", clinDataTmp[,1], "</a>", sep="");
  return(clinDataTmp);
}




####################################################
####################################################
#End Genomically Similar Samples P5
####################################################
####################################################


####################################################
####################################################
#All Findings P6
####################################################
####################################################

allFindingsTable <- function()
{
  #Druggability
  drugData <- filterDruggability();
  
  #First Get Mutations
  tmpMut <- filterMutations();
  tmpMut[,"Abberation"] <- ifelse(tmpMut[,"HGVSp_Short"]!="", paste(tmpMut[,"Hugo_Symbol"], tmpMut[,"HGVSp_Short"], sep=": "), as.character(tmpMut[,"Hugo_Symbol"]));
  tmpMut[,"Type"] <- "Mutation";
  tmpMut[,"Details"] <- paste("Mutation Type: ", tmpMut[,"Variant_Classification"], sep="");
  tmpMut[,"Score"] <- "None";
  tmpMut[,"Trials"] <- "None";
  tmpMut <- merge(tmpMut, drugData, by.x="Hugo_Symbol", by.y="gene_name", all.x=T)
  tmpMut <- tmpMut[,c("Abberation", "Type", "Details", "Score", "Drugs", "Trials")]
  
  #Now Fusions
  tmpFus <- filterFusions_aribba();
  tmpFus[,"Abberation"] <-gsub("--", "-", tmpFus[,"X.fusion_name"]);
  tmpFus[,"Type"] <- "Fusion";
  tmpFus[,"Details"] <- paste("Fusion Type: ",tmpFus[,"Splice_type"], sep="");
  tmpFus[,"Score"] <- "None";
  tmpFus[,"Trials"] <- "None";
  tmpFus <- merge(tmpFus, drugData, by.x="HeadGene", by.y="gene_name", all.x=T)
  tmpFus <- merge(tmpFus, drugData, by.x="TailGene", by.y="gene_name", all.x=T)
  tmpFus[,"Drugs"] <- paste(tmpFus[,"Drugs.x"], tmpFus[,"Drugs.y"], sep=", ")
  tmpFus <- tmpFus[,c("Abberation", "Type", "Details", "Score", "Drugs", "Trials")]

  #Now Copy Number
#  tmpCnv <- filterCNV()
#  tmpCnv[,"Abberation"] <- tmpCnv[,1]
#  tmpCnv[,"Type"] <- ifelse(tmpCnv[,2]>2, "Amplification", "Deletion")
#  tmpCnv[,"Details"] <- paste("Copy Number Value: ",tmpCnv[,2], sep="");
#  tmpCnv[,"Score"] <- "None";
#  tmpCnv[,"Trials"] <- "None";
#  tmpCnv <- merge(tmpCnv, drugData, by.x="Gene", by.y="gene_name", all.x=T)
#  tmpCnv <- tmpCnv[,c("Abberation", "Type", "Details", "Score", "Drugs", "Trials")]
  
  #Now Expression
  tmpExp <- RNASeqAnalysisOut[[1]][[2]]
  tmpExp[,"Abberation"] <-rownames(tmpExp)
  tmpExp[,"Type"] <- c(rep("Outlier-High (mRNA)", 20), rep("Outlier-Low (mRNA)", 20));
  tmpExp[,"Details"] <- paste("Z-Score / FPKM: ",round(tmpExp[,"Z_Score"],2), " / ", tmpExp[,"FPKM"], sep="");
  tmpExp[,"Score"] <- "None";
  tmpExp[,"Trials"] <- "None";
  tmpExp <- merge(tmpExp, drugData, by.x="Abberation", by.y="gene_name", all.x=T)
  tmpExp <- tmpExp[,c("Abberation", "Type", "Details", "Score", "Drugs", "Trials")]

  #Now Pathway
  tmpPath <- RNASeqAnalysisOut[[2]][[2]]
  tmpPathUp <- tmpPath[tmpPath[,"Direction"]=="Up",][1:20,]
  tmpPathDown <- tmpPath[tmpPath[,"Direction"]=="Down",][1:20,]
  tmpPath <- rbind(tmpPathUp, tmpPathDown);
  tmpPath[,"Abberation"] <-gsub("HALLMARK_", "", tmpPath[,"Pathway"]);
  tmpPath[,"Abberation"] <-gsub("_", " ", tmpPath[,"Abberation"]);
  tmpPath[,"Type"] <- c(rep("Pathway Up", 20), rep("Pathway Down", 20));
  tmpPath[,"Details"] <- paste("P-Value: ",as.character(formatC(tmpPath[,"P_VAL"], format = "e", digits = 2)), sep="");
  tmpPath[,"Score"] <- "None";
  tmpPath[,"Drugs"] <- "None";
  tmpPath[,"Trials"] <- "None";
  tmpPath <- tmpPath[,c("Abberation", "Type", "Details", "Score", "Drugs", "Trials")]
  
  #Now Merge Together
#  allFindingsDF <- rbind(tmpMut, tmpFus, tmpCnv, tmpExp, tmpPath);
  allFindingsDF <- rbind(tmpMut, tmpFus, tmpExp, tmpPath);
  allFindingsDF[is.na(allFindingsDF[,"Drugs"]),"Drugs"]<- "None";
  allFindingsDF[allFindingsDF[,"Drugs"]=="NA, NA","Drugs"]<- "None";
  allFindingsDF <- unique(allFindingsDF);
  return(allFindingsDF);
  
}



####################################################
####################################################
#End All Findings P6
####################################################
####################################################



####################################################
####################################################
#Genomic Landscape P7
####################################################
####################################################


####################################################
#Function for CNV View
####################################################
plotCNV <- function(myCnvData=cnvData)
{
  colnames(myCnvData) <- c("Chr", "StartBP", "EndBP", "AbsCNV", "Type");
  myCnvData[,"log2CNV"] <- log2((myCnvData[,"AbsCNV"]+1))
  myCnvData[,"MedianBP"] <- round((myCnvData[,"StartBP"]+myCnvData[,"EndBP"])/2)
  myCnvData <- myCnvData[,c("Chr", "MedianBP", "log2CNV")]
  single.seg <- pcf(data=myCnvData,gamma=12,verbose=FALSE)
  plotGenome(myCnvData, single.seg)
}
####################################################
#End Function for CNV View
####################################################

####################################################
#Function for Circos Plot
####################################################
plotCircos <- function()
{
  data(UCSC.HG38.Human.CytoBandIdeogram); 
  chr.exclude <- NULL;
  cyto.info <- UCSC.HG38.Human.CytoBandIdeogram;
  tracks.inside <- 10;
  tracks.outside <- 0;
  RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside);
  
  rcircos.params <- RCircos.Get.Plot.Parameters();
  rcircos.cyto <- RCircos.Get.Plot.Ideogram();
  rcircos.position <- RCircos.Get.Plot.Positions();
  
  png("tmpRCircos.png")
  RCircos.Set.Plot.Area();  
#  par(mai=c(0.25, 0.25, 0.25, 0.25));
  plot.new();
  plot.window(c(-2.5,2.5), c(-2.5, 2.5));
  RCircos.Chromosome.Ideogram.Plot();

  #Put in Mutations
  mySymbols=as.character(filterMutations()[,1]);
  RCircos.Gene.Label.Data <- chrMap[chrMap[,1]%in%mySymbols,]
  RCircos.Gene.Label.Data <- RCircos.Gene.Label.Data[!grepl("CHR_", RCircos.Gene.Label.Data[,4]), ]
  RCircos.Gene.Label.Data[,4] <- paste("chr", RCircos.Gene.Label.Data[,4], sep="");
  colnames(RCircos.Gene.Label.Data) <- c("hgnc_symbol", "start_position", "end_position", "chromosome_name");
  RCircos.Gene.Label.Data <- RCircos.Gene.Label.Data[,c("chromosome_name", "start_position", "end_position", "hgnc_symbol")]
  RCircos.Gene.Label.Data <- RCircos.Gene.Label.Data[order(RCircos.Gene.Label.Data[,4]),];
  name.col <- 4;
  side <- "in";
  track.num <- 1;
  RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data, track.num, side);
  track.num <- 2;
  RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data, name.col,track.num, side);
  
  #Add Expression
  tmpExp <- data.frame(names(RNASeqAnalysisOut[[1]][[1]]), RNASeqAnalysisOut[[1]][[1]]);
  tmpExp[,2] <- ifelse(tmpExp[,2]>5, 5, ifelse(tmpExp[,2]<(-5), -5, tmpExp[,2]))
  colnames(tmpExp) <- c("GeneName", "Expression")
  RCircos.Heatmap.Data <- merge(chrMap, tmpExp, by.x="HGNC.symbol", by.y="GeneName");
  colnames(RCircos.Heatmap.Data) <- c("GeneName", "chromStart", "chromEnd", "Chromosome", "Expression")
  RCircos.Heatmap.Data <- RCircos.Heatmap.Data[!grepl("CHR_", RCircos.Heatmap.Data[,4]), ]
  RCircos.Heatmap.Data <- RCircos.Heatmap.Data[!grepl("MT", RCircos.Heatmap.Data[,4]), ]
  RCircos.Heatmap.Data[,4] <- paste("chr", RCircos.Heatmap.Data[,4], sep="");
  RCircos.Heatmap.Data <- RCircos.Heatmap.Data[,c("Chromosome", "chromStart", "chromEnd", "GeneName", "Expression")];
  data.col <- 5;
  track.num <- 4;
  side <- "in";
  RCircos.Heatmap.Plot(RCircos.Heatmap.Data, data.col, track.num, side);
  RCircos.Heatmap.Data.High <- RCircos.Heatmap.Data[abs(RCircos.Heatmap.Data[,"Expression"])==5,];
  RCircos.Heatmap.Data.High <- RCircos.Heatmap.Data.High[RCircos.Heatmap.Data.High[,"GeneName"]%in%as.character(cancerGenes[,1]),]
  RCircos.Gene.Connector.Plot(RCircos.Heatmap.Data.High, 5, side);
  RCircos.Gene.Name.Plot(RCircos.Heatmap.Data.High, 4,6, side);
  
  #Add Fusions
  myFus <- filterFusions_aribba()
  RCircos.Link.Data.tmp.h <- chrMap[chrMap[,1]%in%myFus[,"HeadGene"],];
  RCircos.Link.Data.tmp.h <- RCircos.Link.Data.tmp.h[!grepl("CHR_HSCHR1_1_CTG32_1", RCircos.Link.Data.tmp.h[,"Chromosome.scaffold.name"]),]
  RCircos.Link.Data.tmp.h <- RCircos.Link.Data.tmp.h[,c(4,2,3,1)]
  RCircos.Link.Data.tmp.h[,1] <- paste("chr", RCircos.Link.Data.tmp.h[,1], sep="");
  
  RCircos.Link.Data.tmp.t <- chrMap[chrMap[,1]%in%myFus[,"TailGene"],]
  
  RCircos.Link.Data.tmp.t <- RCircos.Link.Data.tmp.t[,c(4,2,3,1)]
  RCircos.Link.Data.tmp.t[,1] <- paste("chr", RCircos.Link.Data.tmp.t[,1], sep="");
  
  RCircos.Link.Data <- data.frame(RCircos.Link.Data.tmp.h[,c(1,2,3)], RCircos.Link.Data.tmp.t[,c(1,2,3)]);
  track.num <- 12;
  RCircos.Link.Plot(RCircos.Link.Data, track.num, TRUE);
  RCircos.Gene.Name.Plot(rbind(RCircos.Link.Data.tmp.h, RCircos.Link.Data.tmp.t), 4,9, inside.pos=50);
  dev.off()
}
plotCircos()


####################################################
#End Function for Circos Plot
####################################################

####################################################
#Function for Network View
####################################################
plotNetwork <- function(numGenes=250)
{
  #Let's build all our nodes
  nodeGenesMut <- as.character(filterMutations()[,1]) #Mutations
  nodeGenesMut <- c(nodeGenesMut, c(filterFusions_aribba()[,"HeadGene"]), c(filterFusions_aribba()[,"TailGene"])) #Fusions
  rnaGenes <-RNASeqAnalysisOut[[1]][[1]]
  rnaGenes <- data.frame(names(rnaGenes), rnaGenes);
  upGenes <- as.character(rnaGenes[order(-rnaGenes[,2]),][1:numGenes,1]);
  downGenes <- as.character(rnaGenes[order(rnaGenes[,2]),][1:numGenes,1]);
  nodeGenes <- c(nodeGenesMut, upGenes, downGenes);
  tmpGeneMania <- geneMania[geneMania[,"Gene_A_EntrezGeneName"]%in%nodeGenes,]
  tmpGeneMania <- tmpGeneMania[tmpGeneMania[,"Gene_B_EntrezGeneName"]%in%nodeGenes,]
  tmpGeneMania <- tmpGeneMania[tmpGeneMania[,"Network_Group_Name"]%in%c("Co-localization", "Genetic Interactions", "Pathway", "Physical Interactions"),]
  cifNetwork <- tmpGeneMania[,c("Gene_A_EntrezGeneName", "Gene_B_EntrezGeneName")];
  colnames(cifNetwork) <- c("from", "to")
  cifNetwork[,"logFC"] <- ifelse(cifNetwork[,"from"]%in%nodeGenesMut, 5, 1)
  cifNetwork[,"miRNA"] <- "No";
  gR<-polishNetwork(cifNetwork)
  tmpNetwork <- igraph.from.graphNEL(gR);
  p <- ggnetwork(tmpNetwork, layout = "fruchtermanreingold", cell.jitter = 0.75);
  nodeSize <- data.frame(table(p[,"vertex.names"]));
  colnames(nodeSize)[2] <- "Frequency"
  p <- merge(p, nodeSize, by.x="vertex.names", by.y="Var1")
  ggplot(p, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(color = "grey50")+
    geom_nodelabel(aes(label = vertex.names, color=Frequency),fontface = "bold")+
    theme_blank()+scale_color_continuous(low="grey", high="red")
}
plotNetwork()
####################################################
#End Function for Network View
####################################################


####################################################
####################################################
#End Genomic Landscape P7
####################################################
####################################################


