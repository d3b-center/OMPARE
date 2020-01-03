######################################
# Universal Function to filter fusions
######################################

# For fusion filtering, use annoFuse gene list instead 
# cancerGenes <- read.delim("data/Reference/CancerGeneList.tsv", stringsAsFactors = F)
cancerGenes <- read.delim("data/Reference/genelistreference.txt", stringsAsFactors = F)
cancerGenes <- subset(cancerGenes, type == "TumorSuppressorGene" | type == "CosmicCensus" | type == "Oncogene")

filterFusions <- function(myFusFile = fusFile, myCancerGenes = cancerGenes, myJunctionReads = 2) {
  
  nm <- grep('star', myFusFile)
  if(length(nm) == 1){
    method = "star"
  } else {
    method = "arriba"
  }
  
  myFusData <- read.delim(myFusFile)
  fusDataFilt <- myFusData
  
  # separate comma separated gene names (arriba)
  # fusDataFilt <- fusDataFilt %>% 
  #   mutate(X.gene1 = strsplit(as.character(X.gene1), ",")) %>% 
  #   unnest(X.gene1) %>%
  #   as.data.frame()
  # fusDataFilt <- fusDataFilt %>% 
  #   mutate(gene2 = strsplit(as.character(gene2), ",")) %>% 
  #   unnest(gene2) %>%
  #   as.data.frame()
  # fusDataFilt$X.gene1 <- gsub('[(].*','',fusDataFilt$X.gene1)
  # fusDataFilt$gene2 <- gsub('[(].*','',fusDataFilt$gene2)
  
  # Format column names
  if(method == "star"){
    fusDataFilt <- fusDataFilt %>%
      mutate(X.fusion_name = gsub('--','_', X.fusion_name),
             HeadGene = gsub('_.*','', X.fusion_name),
             TailGene = gsub('.*_','', X.fusion_name))
  } else {
    fusDataFilt <- fusDataFilt %>%
      mutate(X.fusion_name = paste0(as.character(X.gene1), "_", as.character(gene2)),
             Splice_type = type,
             HeadGene = X.gene1,
             TailGene = gene2)
  }
  
  # Filter by Number of Reads (STAR) or by Confidence (Arriba)
  if(method == "star"){
    fusDataFilt <- fusDataFilt %>% 
      filter(JunctionReads > myJunctionReads)
  } else {
    fusDataFilt <- fusDataFilt %>%
      filter(confidence != "low")
  }
  
  # Filter by Cancer Gene List (Annofuse)
  myCancerGenes <- as.character(myCancerGenes$Gene_Symbol)
  fusDataFilt <- fusDataFilt %>%
    filter(HeadGene %in% myCancerGenes | TailGene %in% myCancerGenes) %>%
    dplyr::select(X.fusion_name, Splice_type, HeadGene, TailGene)
    
  return(fusDataFilt)
}
