####################################
# Tabulate RNA-Seq Pathway Analysis
####################################

# Function to return all results from RNA-Seq Analysis
suppressPackageStartupMessages(library(GSEABase))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RDiseaseXpress))
suppressPackageStartupMessages(library(xlsx))
suppressPackageStartupMessages(library(optparse))

# source code
source('code/helper.R')

option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              help = "Directory e.g. data/PNOC008-04"),
  make_option(c("-o", "--output"), type = "character", 
              help = "output excel file with extension i.e. output.xlsx")
)

# load data
opt <- parse_args(OptionParser(option_list = option_list))
topDir <- opt$input
fname <- opt$output

# Dataset1: GTex Brain
gtexData <- readRDS("data/Reference/GTEx/GTEx_fullExpr_matrix.RDS")
gtexBrain <- getSamples(myStudy = "GTEx")
gtexBrain <- gtexBrain %>%
  filter(subtissue == "Brain")
gtexData <- gtexData[,c("gene_id", gtexBrain$sample_id)]
gtexGeneAnnot <- read.delim("data/Reference/GTEx/gencode.v23.annotation.gi_ti_gs.txt", stringsAsFactors =F)
gtexGeneAnnot <- gtexGeneAnnot[,c("gene_id", "gene_symbol")]

# Format gtex data ensembl ids
gtexData$gene_id <- sapply(as.character(gtexData$gene_id), FUN=remDotStuff)
rownames(gtexData) <- gtexData$gene_id
gtexData <- gtexData[-1]

# PBTA
# clinical data
pbta.hist <- read.delim('~/Projects/OpenPBTA-analysis/data/pbta-histologies.tsv', stringsAsFactors = F)
pbta.hist <- pbta.hist %>%
  filter(experimental_strategy == "RNA-Seq")

# Dataset2: PBTA (full)
pbta.polya <- readRDS('~/Projects/OpenPBTA-analysis/data/pbta-gene-expression-rsem-fpkm.polya.rds')
pbta.stranded <- readRDS('~/Projects/OpenPBTA-analysis/data/pbta-gene-expression-rsem-fpkm.stranded.rds')
pbta.full <- merge(pbta.polya, pbta.stranded, by = 'gene_id')
pbta.full <- pbta.full[grep('PAR_Y', pbta.full$gene_id, invert = TRUE),]
pbta.annot <- pbta.full %>% 
  separate(gene_id, c("gene_id", "gene_symbol"), sep = "\\_", extra = "merge") %>%
  dplyr::select(c(gene_id, gene_symbol))
pbta.full <- pbta.full %>% 
  separate(gene_id, c("gene_id", "gene_symbol"), sep = "\\_", extra = "merge") %>%
  dplyr::select(-c(gene_symbol))
pbta.full$gene_id <- sapply(as.character(pbta.full$gene_id), FUN=remDotStuff)
rownames(pbta.full) <- pbta.full$gene_id
pbta.full <- pbta.full[-1]

# Dataset3: PBTA (HGG)
pbta.hist <- pbta.hist %>%
  filter(disease_type_new == "High-grade glioma")
pbta.hgg <- pbta.full[,colnames(pbta.full) %in% pbta.hist$Kids_First_Biospecimen_ID]

# Cancer Genes
cancerGenes <- read.delim("data/Reference/CancerGeneList.tsv", stringsAsFactors = F)
cancerGenes <- cancerGenes %>%
  filter(Gene_Symbol != "") %>%
  dplyr::select(-Count) %>%
  gather(key = "file", value = "type", -Gene_Symbol) %>%
  mutate(type = file)
geneListRef <- read.delim("data/Reference/genelistreference.txt", stringsAsFactors = F)
geneListRef <- subset(geneListRef, type == "TumorSuppressorGene" | type == "CosmicCensus" | type == "Oncogene")
cancerGenes <- rbind(cancerGenes, geneListRef)
rm(geneListRef)

# Genesets
hallMarkSets <- getGmt("data/Reference/mSigDB/h.all.v6.2.symbols.gmt", collectionType=BroadCollection(), geneIdType= SymbolIdentifier())
hallMarkSets <- geneIds(hallMarkSets)
hallMarkSetsTS <- stack(hallMarkSets)

# Drug Interaction DB
dgidb <- read.delim("data/Reference/DGIdb.txt", stringsAsFactors = F)

# Patient Expression data
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

# function to tabulate DE and Pathway results
runRNASeqAnalysis <- function(expData = NULL, refData = gtexData, refAnnot = gtexGeneAnnot, thresh = 2, comparison) {
  
  # Merge GTEx and Patient data on common genes
  intGenesTmp <- intersect(rownames(refData), rownames(expData))
  mergeDF <- cbind(refData[intGenesTmp,], expData[intGenesTmp,"FPKM"])
  
  # Collapse to unique gene symbols
  # Matrix of GTEx and Patient FPKM data
  refAnnot <- unique(refAnnot[c("gene_id","gene_symbol")])
  refAnnot$gene_id <- sapply(as.character(refAnnot$gene_id), FUN=remDotStuff)
  rownames(refAnnot) <- refAnnot$gene_id
  mergeDF <- cbind(refAnnot[rownames(mergeDF),2], mergeDF)
  mergeDF[,"meanExp"] <- rowMeans(mergeDF[-1]) 
  mergeDF <- mergeDF[order(-mergeDF[,"meanExp"]),]
  mergeDF <- mergeDF[!duplicated(mergeDF[,1]),]
  rownames(mergeDF) <- mergeDF[,1] 
  mergeDF <- mergeDF[-1]
  mergeDF <- mergeDF[-ncol(mergeDF)] 
  colnames(mergeDF)[ncol(mergeDF)] <- "SampleX"
  
  # Calculate Gene Outliers in Patient (top 20 Up and Down)
  getAllOutliers <- function(myMergeDF = mergeDF, getTop = 20, cancerGeneNames = cancerGenes$Gene_Symbol) {
    
    # Filter in Patient: FPKM > 10 
    myMergeDF <- myMergeDF[myMergeDF$SampleX > 10,]
    
    # z-score and return only patient's value
    getZ <- function(x) {
      x <- log2(x+1)
      out <- (x-mean(x))/sd(x)
      return(out[length(out)])
    }
    output <- apply(myMergeDF, FUN=getZ, MARGIN=1)
    
    # full data
    output.df <- data.frame(output, myMergeDF[names(output),"SampleX"])
    colnames(output.df) <- c("Z_Score", "FPKM")
    output.df$DE <- ""
    output.df$DE[which(output.df$Z_Score < (-1*1.5))] <- "Down"
    output.df$DE[which(output.df$Z_Score > 1.5)] <- "Up"
    output.df$Comparison <- comparison
    output.df$CancerGene <- ifelse(rownames(output.df) %in% cancerGeneNames, TRUE, FALSE)
    
    # top 20 only
    outputCanc <- output[intersect(names(output), cancerGeneNames)] # filter to cancer gene list
    outputDown <- sort(outputCanc)[1:getTop] # top 20 down
    outputUp <- sort(outputCanc, T)[1:getTop] # top 20 up
    
    outputUpDF <- data.frame(outputUp, myMergeDF[names(outputUp),"SampleX"])
    colnames(outputUpDF) <- c("Z_Score", "FPKM")
    outputDownDF <- data.frame(outputDown, myMergeDF[names(outputDown),"SampleX"])
    colnames(outputDownDF) <- c("Z_Score", "FPKM")
    
    return(list(expr.genes.z.score = output, 
                diffexpr.top20 = rbind(outputUpDF, outputDownDF),
                output.df = output.df))
  }
  geneAnalysisOut <- getAllOutliers()
  
  # Calculate Pathway Outliers-
  # Currently use Enrichment, but moving forward will use GSVA
  
  # Get Up and Down Genes
  tmpghj <- geneAnalysisOut[[1]]
  upGenes <- names(tmpghj)[which(geneAnalysisOut[[1]] > thresh)]
  downGenes <- names(tmpghj)[which(geneAnalysisOut[[1]] < (-1*thresh))]
  
  # If not enough genes take top 1000 
  if(length(upGenes) < 1000) {
    upGenes <- names(sort(geneAnalysisOut[[1]], decreasing = TRUE))[1:1000]
  }
  if(length(downGenes) < 1000) {
    downGenes <- names(sort(geneAnalysisOut[[1]], decreasing = FALSE))[1:1000]
  }
  
  # Code to run pathway analysis
  runHypGeom <- function(set, genes, n = 20000, universe = NULL) {
    
    if(!is.null(universe)) {
      set <- intersect(set, universe)
    }
    
    # number of white balls
    x <- length(intersect(genes, set))
    ngenes <- as.character(toString(intersect(genes, set)))
    
    # white balls
    m <- length(genes)
    
    # black balls
    n2 <- n-m 
    
    # balls drawn from the urn 
    k <- length(set)
    
    out <- phyper(x-1, m, n2, k, lower.tail=F)
    setSize <- k
    overLap <- x
    numGenes <- m
    
    myRet <- c(setSize, numGenes, overLap, out, ngenes) 
    return(myRet)
  }
  
  # Accessory for functional enrichment
  funcEnrichment <- function(genes, sets, qval=.25, numRet=5, myN=20000, myUniverse=NULL) {
    
    out <- lapply(sets, FUN = runHypGeom, genes = genes, n=myN, universe=myUniverse)
    out <- data.frame(out)
    out <- data.frame(t(out))
    out[,4] <- as.numeric(as.character(out[,4]))
    out$ADJ_P_VAL <- p.adjust(out[,4], method="BH")
    colnames(out)[1:5] <- c("SET_SIZE", "NUM_GENES_INPUT", "OVERLAP", "P_VAL", "GENES")
    out$Pathway <- rownames(out)
    return(out)
  
  }
  
  # up pathways
  upPathways <- funcEnrichment(upGenes, hallMarkSets, qval=1, myN=25000, myUniverse=rownames(mergeDF))
  upPathways <- upPathways %>%
    mutate(Direction = "Up") %>%
    filter(P_VAL < 0.05) %>%
    arrange(P_VAL)
  
  # down pathways
  downPathways <- funcEnrichment(downGenes, hallMarkSets, qval=1, myN=25000, myUniverse=rownames(mergeDF))
  downPathways <- downPathways %>%
    mutate(Direction = "Down") %>%
    filter(P_VAL < 0.05) %>%
    arrange(P_VAL)
  
  # full pathway dataframe
  pathway.df <- rbind(upPathways, downPathways)
  pathway.df <- pathway.df[,c("Pathway", "GENES", "SET_SIZE", "NUM_GENES_INPUT", "OVERLAP", "P_VAL", "ADJ_P_VAL", "Direction")]
  pathway.df$Comparison <- comparison
  
  # expand gene names 
  pathway.df.exp <- pathway.df %>%
    separate_rows(GENES, convert = TRUE) %>%
    dplyr::select(Pathway, GENES)
  
  # now for the full output.df, add drug info and pathway info
  output.df <- geneAnalysisOut$output.df
  output.df$Gene_name <- rownames(output.df)
  dgidb <- dgidb %>%
    mutate(Gene_name = gene_name) %>%
    filter(interaction_types != "" & drug_name != "" & Gene_name != "") %>%
    dplyr::select(Gene_name, drug_name) %>%
    unique() %>%
    group_by(Gene_name) %>% 
    dplyr::summarize(Drugs=paste(drug_name, collapse=", ")) %>%
    as.data.frame()
  output.df <- merge(output.df, dgidb, by = 'Gene_name', all.x = TRUE)
  output.df <- merge(output.df, pathway.df.exp, by.x = 'Gene_name', by.y = 'GENES', all.x = TRUE)
  output.df <- output.df %>%
    group_by(Gene_name) %>%
    mutate(Pathway = toString(Pathway)) %>%
    unique() %>%
    filter(DE != "")
  
  finalOut <- list(pathways = pathway.df, genes = output.df)
  return(finalOut)
}

# output filename
outdir <- file.path(topDir, "Summary")
fname <- file.path(outdir, fname)

# summarize
GTExBrain <- runRNASeqAnalysis(expData = expData, refData = gtexData, refAnnot = gtexGeneAnnot, comparison = paste0("GTExBrain_", ncol(gtexData)))
PBTA_HGG <- runRNASeqAnalysis(expData = expData, refData = pbta.hgg, refAnnot = pbta.annot, comparison = paste0("PBTA_HGG_", ncol(pbta.hgg)))
PBTA_All <- runRNASeqAnalysis(expData = expData, refData = pbta.full, refAnnot = pbta.annot, comparison = paste0("PBTA_All_", ncol(pbta.full)))

# create one file per output
pathway.df <- rbind(GTExBrain$pathways, PBTA_HGG$pathways, PBTA_All$pathways)
pathway.df <- pathway.df %>%
  group_by(Pathway, Direction) %>%
  mutate(Freq = n()) %>%
  as.data.frame()
pathway.df.up <- pathway.df %>%
  filter(Direction == "Up") %>%
  as.data.frame()
pathway.df.down <- pathway.df %>%
  filter(Direction == "Down") %>%
  as.data.frame()

genes.df <- rbind(GTExBrain$genes, PBTA_HGG$genes, PBTA_All$genes)
genes.df <- genes.df %>%
  group_by(Gene_name, DE) %>%
  mutate(Freq = n(), Drugs = replace_na(Drugs, "NA")) %>% 
  as.data.frame()
genes.df.up <- genes.df %>%
  filter(DE == "Up") %>%
  as.data.frame()
genes.df.down <- genes.df %>%
  filter(DE == "Down") %>%
  as.data.frame()

# write out to excel
write.xlsx(x = pathway.df.up, file = fname, sheetName = "Pathways_Up", row.names = F)
write.xlsx(x = genes.df.up, file = fname, sheetName = "DE_Genes_Up", row.names = F, append = TRUE)
write.xlsx(x = pathway.df.down, file = fname, sheetName = "Pathways_Down", row.names = F, append = TRUE)
write.xlsx(x = genes.df.down, file = fname, sheetName = "DE_Genes_Down", row.names = F, append = TRUE)
