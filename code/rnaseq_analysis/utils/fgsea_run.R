# function for fgsea enrichment
fgsea_run <- function(pathways, ranks, comparison){
  set.seed(42)
  fgseaRes <- fgsea(pathways, ranks, minSize = 15, maxSize = 1500, eps = 0.0)
  
  # significant pathways
  fgseaRes <- fgseaRes %>%
    filter(padj < 0.05) %>%
    mutate(direction = ifelse(ES > 0, "up", "down")) %>%
    arrange(padj) %>%
    mutate(genes = sapply(leadingEdge, paste, collapse=",")) %>%
    mutate(comparison = comparison)  %>%
    dplyr::select(-c(log2err, leadingEdge))
  
  return(fgseaRes)
}