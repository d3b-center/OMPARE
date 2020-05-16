# Author: Komal S. Rathi
# Date: 04/25/2020
# Function: Up/Down pathways for each PBTA sample, compare to rest of PBTA (1) and GTEx (2)
# do this once and read in for tabulate pathways (Page 8 of report)

# Function to return all results from RNA-Seq Analysis
suppressPackageStartupMessages(library(GSEABase))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(xlsx))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(doMC))
registerDoMC(cores = 4)

# Dataset1: GTex Brain
# clinical
if(!file.exists('data/Reference/GTEx/GTEx_clinical.RDS')){
  gtexBrain <- getSamples(myStudy = "GTEx")
} else {
  gtexBrain <- readRDS('data/Reference/GTEx/GTEx_clinical.RDS')
}
gtexBrain <- gtexBrain %>%
  filter(subtissue == "Brain")

# gencode annotation
gtexGeneAnnot <- read.delim("data/Reference/GTEx/gencode.v23.annotation.gi_ti_gs.txt", stringsAsFactors =F)

# expression
gtexData <- readRDS("data/Reference/GTEx/GTEx_fullExpr_matrix.RDS")
gtexData <- gtexData[,c("gene_id", gtexBrain$sample_id)]
gtexData <- gtexGeneAnnot %>%
  dplyr::select(c(gene_id, gene_symbol)) %>%
  inner_join(gtexData, by = 'gene_id') %>%
  mutate(means = rowMeans(dplyr::select(.,-gene_id, -gene_symbol))) %>% # take rowMeans
  arrange(desc(means)) %>% # arrange decreasing by means
  distinct(gene_symbol, .keep_all = TRUE) %>% # keep the ones with greatest mean value. If ties occur, keep the first occurencce
  dplyr::select(-c(gene_id, means)) %>%
  column_to_rownames('gene_symbol') 

# Dataset2: PBTA (full n = 970)
# clinical
pbta.hist <- read.delim('data/Reference/PBTA/pbta-histologies.tsv', stringsAsFactors = F)
pbta.hist <- pbta.hist %>%
  filter(experimental_strategy == "RNA-Seq",
         RNA_library == "stranded")

# expression
pbta.stranded <- readRDS('data/Reference/PBTA/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds')
pbta.full <- pbta.stranded

# Dataset3: PBTA (HGG n = 96)
pbta.hist <- pbta.hist %>%
  filter(integrated_diagnosis == "High-grade glioma")
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
dgidb <- dgidb %>%
  mutate(Gene_name = gene_name) %>%
  filter(interaction_types != "" & drug_name != "" & Gene_name != "") %>%
  dplyr::select(Gene_name, drug_name) %>%
  unique() %>%
  group_by(Gene_name) %>%
  dplyr::summarize(Drugs=paste(drug_name, collapse=", ")) %>%
  as.data.frame()

# function to tabulate DE and Pathway results
runRNASeqAnalysis <- function(expData = NULL, refData = gtexData, thresh = 2, comparison) {
  
  smp <- unique(as.character(expData$Sample))
  
  # done for PBTA
  expData <- expData %>%
    remove_rownames() %>%
    column_to_rownames('Gene') %>%
    dplyr::select(-c(Sample))
  
  # now remove the current PBTA SOI from refData
  refData <- refData[,grep(smp, colnames(refData), invert = T, value = T)]
  
  # Merge GTEx and Patient data on common genes
  intGenesTmp <- intersect(rownames(refData), rownames(expData))
  mergeDF <- cbind(refData[intGenesTmp,], expData[intGenesTmp, "FPKM"])
  
  # Collapse to unique gene symbols
  # Matrix of reference data and Patient FPKM data
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
    filter(ADJ_P_VAL < 0.05) %>%
    arrange(ADJ_P_VAL)
  
  # down pathways
  downPathways <- funcEnrichment(downGenes, hallMarkSets, qval=1, myN=25000, myUniverse=rownames(mergeDF))
  downPathways <- downPathways %>%
    mutate(Direction = "Down") %>%
    filter(ADJ_P_VAL < 0.05) %>%
    arrange(ADJ_P_VAL)
  
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

res <- melt(as.matrix(pbta.full), value.name = "FPKM", varnames = c("Gene", "Sample"))
system('mkdir -p data/Reference/GSEA')
if(!file.exists('data/Reference/GSEA/PBTA_vs_GTExBrain.RDS')){
  GTExBrain <-  plyr::dlply(res, .variables = "Sample", .fun = function(x) runRNASeqAnalysis(expData = x, refData = gtexData, comparison = paste0("GTExBrain_", ncol(gtexData))), .parallel = TRUE)
  saveRDS(GTExBrain, file = 'data/Reference/GSEA/PBTA_vs_GTExBrain.RDS')
}

if(!file.exists('data/Reference/GSEA/PBTA_vs_PBTA.RDS')){
  PBTA_All <- plyr::dlply(res, .variables = "Sample", .fun = function(x) runRNASeqAnalysis(expData = x, refData = pbta.full, comparison = paste0("PBTA_All_", ncol(pbta.full))), .parallel = TRUE)
  saveRDS(PBTA_All, file = 'data/Reference/GSEA/PBTA_vs_PBTA.RDS')
}

if(!file.exists('data/Reference/GSEA/PBTA_vs_PBTAHGG.RDS')){
  PBTA_HGG <- plyr::dlply(res, .variables = "Sample", .fun = function(x) runRNASeqAnalysis(expData = x, refData = pbta.hgg, comparison = paste0("PBTA_HGG_", ncol(pbta.hgg))), .parallel = TRUE)
  saveRDS(PBTA_HGG, file = 'data/Reference/GSEA/PBTA_vs_PBTAHGG.RDS')
}
