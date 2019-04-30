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


source("code/RNASeqAnalysis.R")
source("code/dgidbParse.R")

source("code/copyNumberConvert.R")
source("code/patientSampleInfo.R")
source("code/pubTheme.R");
#############################
#-End Load packages & source
#############################


#############################
#-Read Data
#############################
#Mutation Data
mutData <- read.delim("../data/CBTTC-HGG/MutationsMAF/5c3eace5-950a-4a05-81ed-5c04b4a0a367.strelka.vep.maf", skip=1);

#Copy Number Data
cnvData <- read.delim("../data/CBTTC-HGG/CNV/5c3eace5-950a-4a05-81ed-5c04b4a0a367.CNVs", header=F);
cnvGenes <- createCopyNumber()

#Expression Data
expData <- read.delim("../data/CBTTC-HGG/ExpressionGene/7316-37_564442.genes.results")

#Fusion Data
fusData <- read.delim("../data/CBTTC-HGG/Fusions/7316-37.local.transcript.converted.pe.star-fusion.fusion_candidates.final")

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
mapping <- read.delim("../data/Reference/mappingFile.txt", header=F, stringsAsFactors=F);
clinData <- read.delim("../data/Reference/study_view_clinical_data.txt", stringsAsFactors=F);
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
filterFusions <- function(myFusData=fusData, myCancerGenes=cancerGenes, myJunctionReads=5)
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
                                    ifelse(is.na(myTableFus[,"ind.x"]), "", paste(", Pathway: ", myTableFus[,"ind"], ",", sep="")),
                                    ifelse(is.na(myTableFus[,"ind.y"]), "", paste("", myTableFus[,"ind"], sep="")), sep="")
  myTableFus <- myTableFus[,c("Abberation", "Type", "Details", "Score", "Drugs", "Trials", "SupportEv")]
  
  
  myTable <- rbind(myTableAmp, myTableDel, myTableMut, myTableFus)
  colnames(myTable)[ncol(myTable)] <- "Supporting Evidence";
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
  p <- p+xlab("Gene Symbol")+scale_fill_manual(values = c("red", "blue"));
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
  pathData <- rbind(pathDataUp, pathDataDown)
  
  p <- ggplot(pathData, aes(factor(Pathway), y=(-1)*log10(P_VAL), fill=Direction))+geom_bar(stat="identity")+coord_flip()+theme_bw();
  p <- p+xlab("Pathway Name")+scale_fill_manual(values = c("red", "blue"))+ylab("-log10 P-Value");
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
#Immune/Tumor Signatures P3
####################################################
####################################################


ImmuneProfile <- function()
{
#  resAll[,"max"] <- apply(resAll[3:ncol(res)], FUN=max, MARGIN=1)
#  resAll <- resAll[order(-resAll[,"max"]),]
#  resAll <- resAll[!duplicated(resAll[,2]),]
#  rownames(resAll) <- resAll[,2];
#  resAll <- resAll[-1:-2];
#  resAll <- resAll[-ncol(resAll)]; #Remove max

#raw.scores = rawEnrichmentAnalysis(as.matrix(resAll),
#                                   xCell.data$signaturessca,
#                                   xCell.data$genes)

#  raw.scores[,"CellType"] <- rownames(raw.scores)

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
#Genomically Similar Samples P4
####################################################
####################################################

clinData <- merge(clinData, mapping, by.x="Sample.ID", by.y="V1")
clinDataOrig <- clinData;
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
res <- cbind(res, expData[,"FPKM"]);
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
  p <- p+scale_color_manual(breaks = as.character(tsneData[,"Cancer.Type"]), values=distinctColorPalette(length(unique(tmpCol))))
  p <- p+theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
               legend.title=element_text(size=12), 
               legend.text=element_text(size=12))+labs(colour="Cancer Type")+ guides(shape=FALSE, size=FALSE)
  p;  
}

allCor <- cor(resTmp[ncol(resTmp)], resTmp)
allCor <- data.frame(t(allCor))
allCor[,"samps"] <- rownames(allCor);
allCor <- allCor[order(-allCor[,1]),]
allCor <- allCor[!grepl("Pat", rownames(allCor)),]

getKMPlot <- function(numNeighbors=200)
{
  mySamps <- allCor[1:numNeighbors,"samps"]
  survData[,"group"] <- survData[,"samps"]%in%mySamps
  survData[,"group"] <- ifelse(survData[,"group"], "Cluster With Patient", "Cluster away from Patient")
  fit <- survfit(Surv(OS_Time, OS_Event) ~ group, data = survData)
  ggsurvplot(
    fit,                     # survfit object with calculated statistics.
    data = survData,  # data used to fit survival curves. 
    risk.table = TRUE,       # show risk table.
    pval = TRUE,             # show p-value of log-rank test.
    conf.int = TRUE,         # show confidence intervals for 
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


getSimilarPatients <- function(numNeighbors=10)
{
  mySamps <- allCor[1:numNeighbors,"samps"]
  clinDataTmp <- clinDataOrig[clinDataOrig[,"V2"]%in%mySamps,]
  clinDataTmp <- clinDataTmp[,c("Sample.ID", "Cancer.Type", "AGE", "ETHNICITY", "RACE", "TUMOR_SITE")]
  clinDataTmp[,"Report"] <- paste("<a href='www.google.com'>", clinDataTmp[,1], "</a>", sep="");
  return(clinDataTmp);
}




####################################################
####################################################
#End Genomically Similar Samples P4
####################################################
####################################################


####################################################
####################################################
#All Findings P5
####################################################
####################################################

allFindingsTable <- function()
{
  #Druggability
  drugData <- filterDruggability();
  
  #First Get Mutations
  tmpMut <- filterMutations();
  tmpMut[,"Abberation"] <- paste(tmpMut[,"Hugo_Symbol"], tmpMut[,"HGVSp_Short"], sep=": ");
  tmpMut[,"Type"] <- "Mutation";
  tmpMut[,"Details"] <- paste("Mutation Type: ", tmpMut[,"Variant_Classification"], sep="");
  tmpMut[,"Score"] <- "None";
  tmpMut[,"Trials"] <- "None";
  tmpMut <- merge(tmpMut, drugData, by.x="Hugo_Symbol", by.y="gene_name", all.x=T)
  tmpMut <- tmpMut[,c("Abberation", "Type", "Details", "Score", "Drugs", "Trials")]
  
  #Now Fusions
  tmpFus <- filterFusions();
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
  tmpCnv <- filterCNV()
  tmpCnv[,"Abberation"] <- tmpCnv[,1]
  tmpCnv[,"Type"] <- ifelse(tmpCnv[,2]>2, "Amplification", "Deletion")
  tmpCnv[,"Details"] <- paste("Copy Number Value: ",tmpCnv[,2], sep="");
  tmpCnv[,"Score"] <- "None";
  tmpCnv[,"Trials"] <- "None";
  tmpCnv <- merge(tmpCnv, drugData, by.x="Gene", by.y="gene_name", all.x=T)
  tmpCnv <- tmpCnv[,c("Abberation", "Type", "Details", "Score", "Drugs", "Trials")]
  
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
  allFindingsDF <- rbind(tmpMut, tmpFus, tmpCnv, tmpExp, tmpPath);
  allFindingsDF[is.na(allFindingsDF[,"Drugs"]),"Drugs"]<- "None";
  allFindingsDF[allFindingsDF[,"Drugs"]=="NA, NA","Drugs"]<- "None";
  return(allFindingsDF);
  
}



####################################################
####################################################
#End All Findings P5
####################################################
####################################################



####################################################
####################################################
#Genomic Landscape P6
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
  myFus <- filterFusions()
  RCircos.Link.Data.tmp.h <- chrMap[chrMap[,1]%in%myFus[,"HeadGene"],];
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
  nodeGenesMut <- c(nodeGenesMut, c(filterFusions()[,"HeadGene"]), c(filterFusions()[,"TailGene"])) #Fusions
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
#End Genomic Landscape P6
####################################################
####################################################


