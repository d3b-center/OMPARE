# filter cnvs by oncogenes/tsgs as well as copy number status cutoff

filter_cnv <- function(myCNVData = cnvGenes, myCancerGenes = cancer_genes) {

  # tumor suppressor genes
  myTSGenes <- myCancerGenes %>%
    filter(type %in% c("TumorSuppressorGene", "Is.Tumor.Suppressor.Gene")) %>%
    .$Gene_Symbol %>%
    unique()
  
  # oncogenes
  myOncogenes <- myCancerGenes %>%
    filter(type %in% c("Oncogene", "OncoKB.Annotated", "Is.Oncogene")) %>%
    .$Gene_Symbol %>%
    unique()
  myOncogenes <- setdiff(myOncogenes, myTSGenes)
  
  # gain in oncogenes
  cnvDataFiltUp <- myCNVData %>%
    filter(status == "gain" & hgnc_symbol %in% myOncogenes)
  
  # loss in tsgs
  cnvDataFiltDown <- myCNVData %>%
    filter(status == "loss" & hgnc_symbol %in% myTSGenes)
  
  # combine
  cnvDataFilt <- rbind(cnvDataFiltUp, cnvDataFiltDown)
  
  return(cnvDataFilt)
}
