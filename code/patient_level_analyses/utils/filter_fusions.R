# filter fusions

filter_fusions <- function(myFusFile = fusFile, myCancerGenes = cancerGenes, myJunctionReads = 2) {
  
  nm <- grep('star', myFusFile)
  if(length(nm) == 1){
    method = "star"
  } else {
    method = "arriba"
  }
  
  myFusData <- read.delim(myFusFile)
  fusDataFilt <- myFusData
  
  # format column names
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
  
  # filter by number of reads (star) or confidence (arriba)
  if(method == "star"){
    fusDataFilt <- fusDataFilt %>% 
      filter(JunctionReads > myJunctionReads)
  } else {
    fusDataFilt <- fusDataFilt %>%
      filter(confidence != "low")
  }
  
  # filter by cancer genes (from annoFuse)
  myCancerGenes <- as.character(myCancerGenes$Gene_Symbol)
  fusDataFilt <- fusDataFilt %>%
    filter(HeadGene %in% myCancerGenes | TailGene %in% myCancerGenes) %>%
    dplyr::select(X.fusion_name, Splice_type, HeadGene, TailGene)
    
  return(fusDataFilt)
}
