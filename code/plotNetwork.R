###########################
# Function for Network View
###########################

plotNetwork <- function(numGenes = 250) {
  nodeGenesMut <- ''
  # Let's build all our nodes
  if(exists('mutData')){
    nodeGenesMut <- as.character(filterMutations()[,1]) #Mutations
  }
  if(exists('fusData')){
    nodeGenesMut <- c(nodeGenesMut, c(fusData[,"HeadGene"]), c(fusData[,"TailGene"])) #Fusions
  }
  rnaGenes <-RNASeqAnalysisOut[[1]][[1]]
  rnaGenes <- data.frame(names(rnaGenes), rnaGenes)
  upGenes <- as.character(rnaGenes[order(-rnaGenes[,2]),][1:numGenes,1])
  downGenes <- as.character(rnaGenes[order(rnaGenes[,2]),][1:numGenes,1])
  nodeGenes <- c(nodeGenesMut, upGenes, downGenes)
  tmpGeneMania <- geneMania[geneMania[,"Gene_A_EntrezGeneName"]%in%nodeGenes,]
  tmpGeneMania <- tmpGeneMania[tmpGeneMania[,"Gene_B_EntrezGeneName"]%in%nodeGenes,]
  tmpGeneMania <- tmpGeneMania[tmpGeneMania[,"Network_Group_Name"]%in%c("Co-localization", "Genetic Interactions", "Pathway", "Physical Interactions"),]
  cifNetwork <- tmpGeneMania[,c("Gene_A_EntrezGeneName", "Gene_B_EntrezGeneName")]
  colnames(cifNetwork) <- c("from", "to")
  cifNetwork[,"logFC"] <- ifelse(cifNetwork[,"from"]%in%nodeGenesMut, 5, 1)
  cifNetwork[,"miRNA"] <- "No"
  gR<-polishNetwork(cifNetwork)
  tmpNetwork <- igraph.from.graphNEL(gR)
  p <- ggnetwork(tmpNetwork, layout = "fruchtermanreingold", cell.jitter = 0.75)
  nodeSize <- data.frame(table(p[,"vertex.names"]))
  colnames(nodeSize)[2] <- "Frequency"
  p <- merge(p, nodeSize, by.x="vertex.names", by.y="Var1")
  ggplot(p, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(color = "grey50")+
    geom_nodelabel(aes(label = vertex.names, color=Frequency),fontface = "bold")+
    theme_blank()+scale_color_continuous(low="grey", high="red")
}