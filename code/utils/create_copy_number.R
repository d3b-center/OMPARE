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
    unique()

  return(output)
}










