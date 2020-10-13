network_plot <- function(numGenes = 250, geneMania) {
  
  # build network nodes
  nodeGenesMut <- ''
  
  # filtered mutations
  if(exists('mutDataFilt')){
    nodeGenesMut <- unique(mutDataFilt$Hugo_Symbol) 
  }
  # filtered fusions
  if(exists('fusData')){
    nodeGenesMut <- c(nodeGenesMut, c(fusData$HeadGene), c(fusData$TailGene)) 
  }
  # expression z-scores
  if(exists('expData')){
    rnaGenes <- rnaseq_analysis_output$expr.genes.z.score
    rnaGenes <- data.frame(gene = names(rnaGenes), z.score = rnaGenes)
    upGenes <- rnaGenes %>% arrange(desc(z.score)) %>% slice_head(n = numGenes) %>% .$gene
    downGenes <- rnaGenes %>% arrange(z.score) %>% slice_head(n = numGenes)  %>% .$gene
    nodeGenes <- c(nodeGenesMut, upGenes, downGenes)
  } else {
    nodeGenes <- c(nodeGenesMut)
  }
  
  # subset network objects
  nodeGenes <- unique(nodeGenes)
  cifNetwork <- geneMania %>%
    filter(Gene_A_EntrezGeneName %in% nodeGenes,
           Gene_B_EntrezGeneName %in% nodeGenes,
           Network_Group_Name %in% c("Co-localization", "Genetic Interactions", "Pathway", "Physical Interactions")) %>%
    mutate(from = Gene_A_EntrezGeneName, to = Gene_B_EntrezGeneName,
           logFC = ifelse(from %in% nodeGenesMut, 5, 1),
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