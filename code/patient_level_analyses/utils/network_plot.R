network_plot <- function(transcriptome_drug_rec_output, geneMania) {
  
  # filter to GTEx only
  transcriptome_drug_rec_output <- transcriptome_drug_rec_output %>%
    filter(Comparison == "GTExBrain_1152")
  
  # filtered mutations
  if(exists('mutDataFilt')){
    nodeGenesMut <- unique(mutDataFilt$Hugo_Symbol) 
  } else {
    nodeGenesMut <- ""
  }
  
  # filtered fusions
  if(exists('fusData')){
    fusion_head <- fusData %>%
      mutate(HeadGene = strsplit(as.character(HeadGene), ",")) %>% 
      unnest(HeadGene) %>%
      mutate(HeadGene = gsub("[(].*", "", HeadGene)) %>%
      .$HeadGene %>% unique()
    fusion_tail <- fusData %>%
      mutate(TailGene = strsplit(as.character(TailGene), ",")) %>% 
      unnest(TailGene) %>%
      mutate(TailGene = gsub("[(].*", "", TailGene)) %>%
      .$TailGene %>% unique()
    nodeGenesFus <- c(fusion_head, fusion_tail) 
  } else {
    nodeGenesFus <- ""
  }
  
  # expression logfc
  if(exists('expData')){
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
