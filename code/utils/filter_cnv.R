# filter cnvs by oncogenes/tsgs as well as copy number status cutoff

filter_cnv <- function(myCNVData, myCancerGenes = cancer_genes) {

  # tumor suppressor genes
  myTSGenes <- myCancerGenes %>%
    filter(grepl(pattern = "TumorSuppressorGene|Is.Tumor.Suppressor.Gene", type)) %>%
    .$Gene_Symbol %>%
    unique()
  
  # oncogenes
  myOncogenes <- myCancerGenes %>%
    filter(grepl(pattern = "Oncogene|OncoKB.Annotated|Is.Oncogene", type)) %>%
    .$Gene_Symbol %>%
    unique()
  myOncogenes <- setdiff(myOncogenes, myTSGenes)
  
  # gain in oncogenes
  cnvDataFiltUp <- myCNVData %>%
    filter(status %in% c("Gain", "Amplification") & hgnc_symbol %in% myOncogenes)
  
  # loss in tsgs
  cnvDataFiltDown <- myCNVData %>%
    filter(status %in% c("Loss", "Complete Loss") & hgnc_symbol %in% myTSGenes)
  
  # combine
  cnvDataFilt <- rbind(cnvDataFiltUp, cnvDataFiltDown)
  
  return(cnvDataFilt)
}
