######################################
# Universal Function to filter fusions
######################################

cancerGenes <- read.delim("data/Reference/CancerGeneList.tsv", stringsAsFactors = F)

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
  fusDataFilt <- fusDataFilt %>% 
    mutate(X.gene1 = strsplit(as.character(X.gene1), ",")) %>% 
    unnest(X.gene1) %>%
    as.data.frame()
  fusDataFilt <- fusDataFilt %>% 
    mutate(gene2 = strsplit(as.character(gene2), ",")) %>% 
    unnest(gene2) %>%
    as.data.frame()
  fusDataFilt$X.gene1 <- gsub('[(].*','',fusDataFilt$X.gene1)
  fusDataFilt$gene2 <- gsub('[(].*','',fusDataFilt$gene2)
  
  if(method == "star"){
    fusDataFilt[,'X.fusion_name'] <- gsub('--','_',fusDataFilt[,'X.fusion_name'])
    fusDataFilt[,"HeadGene"] <- gsub('_.*','',fusDataFilt[,'X.fusion_name'])
    fusDataFilt[,"TailGene"] <- gsub('.*_','',fusDataFilt[,'X.fusion_name'])
  } else {
    fusDataFilt[,"X.fusion_name"] <- paste0(as.character(fusDataFilt[,"X.gene1"]), "_", as.character(fusDataFilt[,"gene2"]))
    fusDataFilt[,"Splice_type"] <- fusDataFilt[,"type"]
    fusDataFilt[,"HeadGene"] <- as.character(fusDataFilt[,"X.gene1"])
    fusDataFilt[,"TailGene"] <- as.character(fusDataFilt[,"gene2"])
  }
  
  # Filter by Number of Reads
  if(method == "star"){
    fusDataFilt <- fusDataFilt[fusDataFilt[,"JunctionReads"] > myJunctionReads,]
  } else {
    fusDataFilt <- fusDataFilt[fusDataFilt[,"confidence"] != "low",]
  }
  
  # Filter by Cancer Gene Census
  myCancerGenes <- as.character(myCancerGenes[,1])
  fusDataFilt <- fusDataFilt[fusDataFilt[,"HeadGene"] %in% myCancerGenes | fusDataFilt[,"TailGene"] %in% myCancerGenes,]
  
  # Subset to required columns only
  fusDataFilt <- fusDataFilt[,c('X.fusion_name','Splice_type','HeadGene','TailGene')]
  return(fusDataFilt)
}
