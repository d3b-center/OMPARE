suppressPackageStartupMessages({
  library(GeneNetworkBuilder)
  library(igraph)
  library(ggnetwork)
})

network_plot <- function(transcriptome_drug_rec_output, geneMania, filtered_mutations, filtered_fusions, rnaseq_analysis_output) {
  
  # filter to GTEx only
  transcriptome_drug_rec_output <- transcriptome_drug_rec_output %>%
    filter(grepl("Normal", Comparison))
  
  # filtered mutations
  if(!is.null(filtered_mutations)){
    nodeGenesMut <- unique(filtered_mutations$Hugo_Symbol) 
  } else {
    nodeGenesMut <- ""
  }
  
  # filtered fusions
  if(!is.null(filtered_fusions)){
    fusion_head <- filtered_fusions %>%
      mutate(gene1 = strsplit(as.character(gene1), ",")) %>% 
      unnest(gene1) %>%
      mutate(gene1 = gsub("[(].*", "", gene1)) %>%
      .$gene1 %>% unique()
    fusion_tail <- filtered_fusions %>%
      mutate(gene2 = strsplit(as.character(gene2), ",")) %>% 
      unnest(gene2) %>%
      mutate(gene2 = gsub("[(].*", "", gene2)) %>%
      .$gene2 %>% unique()
    nodeGenesFus <- c(fusion_head, fusion_tail) 
  } else {
    nodeGenesFus <- ""
  }
  
  # expression logfc
  if(!is.null(rnaseq_analysis_output)){
    rnaGenes <- rnaseq_analysis_output$expr.genes.logfc
    rnaGenes <- data.frame(gene = names(rnaGenes), logfc = rnaGenes)
    rnaGenes <- rnaGenes %>%
      filter(gene %in% transcriptome_drug_rec_output$Gene)
    nodeGenesExp <- unique(rnaGenes$gene)
  } else {
    nodeGenesExp <- ""
  }
  
  # subset network objects
  nodeGenes <- unique(c(nodeGenesMut, nodeGenesFus, nodeGenesExp))
  cifNetwork <- geneMania %>%
    filter(Gene_A_EntrezGeneName %in% nodeGenes,
           Network_Group_Name %in% c("Co-localization", "Genetic Interactions", "Pathway", "Physical Interactions")) %>%
    mutate(from = Gene_A_EntrezGeneName, 
           to = Gene_B_EntrezGeneName) %>%
    left_join(rnaGenes, by = c("from" = "gene")) %>%
    mutate(logFC = ifelse(from %in% nodeGenesMut, 5, logfc),
           logFC = ifelse(is.na(logFC), 1, logFC),
           miRNA = "No") %>%
    dplyr::select(from, to, logFC, miRNA)
  
  # plot network    
  set.seed(100)
  gR <- polishNetwork(cifNetwork)
  tmpNetwork <- igraph.from.graphNEL(gR)
  tmpNetwork <- ggnetwork(tmpNetwork, layout = "fruchtermanreingold", cell.jitter = 0.75)
  nodeSize <- data.frame(table(tmpNetwork$vertex.names))
  tmpNetwork <- tmpNetwork %>%
    inner_join(nodeSize, by = c("vertex.names" = "Var1")) %>%
    dplyr::rename("Frequency" = "Freq") 
  p <- ggplot(tmpNetwork, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(color = "grey50")+
    geom_nodelabel(aes(label = vertex.names, color = Frequency), fontface = "bold")+
    theme_blank() + 
    scale_color_continuous(low = "grey", high = "red")
  return(p)
}
