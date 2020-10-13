#######################################
# Copy number by gene instead of region
#######################################

suppressPackageStartupMessages(library(scales))

create_copy_number <- function(cnvData, ploidy = NULL){
  
  # chromosome map
  chr_map <- chr_map %>%
    filter(hgnc_symbol != "")
  
  # map to cnvData
  output <- cnvData %>%
    mutate(WilcoxonRankSumTestPvalue = as.numeric(scales::scientific(WilcoxonRankSumTestPvalue, digits = 3))) %>%
    rowwise() %>%
    inner_join(chr_map, by = c("chr" = "chromosome")) %>%
    filter(gene_start > start, 
           gene_end < end) %>%
    dplyr::select(hgnc_symbol, copy.number, status, WilcoxonRankSumTestPvalue) %>%
    unique()
  
  # for missing genes
  # if ploidy is not supplied explicitly, use value corresponding to neutral status
  if(is.null(ploidy)){
    ploidy <- output %>%
      filter(status  == "neutral") %>%
      .$copy.number %>%
      min()
  }

  # if there is no CN entry for neutral status, then do this
  if(length(ploidy) == 0){
    ploidy <- output %>%
      filter(status == "gain") %>%
      .$copy.number %>%
      min-1
  }
  missingGenes <- setdiff(chr_map$hgnc_symbol, output$hgnc_symbol)
  missingGenes <- data.frame(missingGenes, ploidy, 'neutral', 1)
  colnames(missingGenes) <- colnames(output)
  output <- unique(rbind(output, missingGenes))

  return(output)
}










