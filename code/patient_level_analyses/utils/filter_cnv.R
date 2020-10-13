# filter cnvs by oncogenes/tsgs as well as copy number status cutoff

filter_cnv <- function(myCNVData = cnvGenes, myCancerGenes = cancerGenes) {

  myTSGenes <- myCancerGenes %>%
    filter(type == "TumorSuppressorGene") %>%
    .$Gene_Symbol %>%
    unique()
  myOncogenes <- myCancerGenes %>%
    filter(type %in% c("Oncogene", "OncoKB")) %>%
    .$Gene_Symbol %>%
    unique()
  myOncogenes <- setdiff(myOncogenes, myTSGenes)
  cnvDataFiltUp <- myCNVData %>%
    filter(status == "gain" & hgnc_symbol %in% myOncogenes)
  cnvDataFiltDown <- myCNVData %>%
    filter(status == "loss" & hgnc_symbol %in% myTSGenes)
  cnvDataFilt <- rbind(cnvDataFiltUp, cnvDataFiltDown)
  
  return(cnvDataFilt)
}
